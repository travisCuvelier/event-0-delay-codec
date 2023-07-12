classdef sortedAdaptiveCutoffPMF64 < handle
         

    properties (SetAccess = public, GetAccess = public)
        total_iterations
        counts %actual number of symbols
        maxCount 
        nSymbols
        nSymbolsPerDim
        nDims
        perm
        iperm
    end
    
    properties(Constant)
        wordSize = 32; %cdf cannot exceed 2^wordSize-1
        precision = 64;     
    end
    
    methods(Static=true)
        
        %big endian binary expansion of z (assumed column) given bytewidth.
        function ze = binaryExpand(z,bytewidth)
            ze = zeros(numel(z),bytewidth);
            
            for idx = 1:bytewidth
                ze(:,bytewidth+1-idx) = mod(z,2);
                z = z - ze(:,bytewidth+1-idx);
                z = z/2;
            end
        end
        
        %this is like bi2de(,'left-msb') in the communications toolbox.
        function d = binary2decimal(b)
            
            [~,nc] = size(b);
            
            d = sum((2.^((nc-1):-1:0)).*b,2);
            
        end %the rate of the code in bits/channel use
        
    end
    
    
    methods
        %constructor
        %precondition (we do not check this)
        %prod(1+biggestSymbolPerDimension) <= 2^64-1. If you do not do
        %this, stuff will go poorly for you. 
        function obj = sortedAdaptiveCutoffPMF64(counts,biggestSymbolPerDimension)

            if(~isrow(counts))
                counts = counts.';
                if(~isrow(counts))
                    error('counts should be a one dimensional array')
                end
            end

            obj.nSymbols = uint64(numel(counts));
            
            if(obj.nSymbols~=prod(uint64(biggestSymbolPerDimension)+uint64(1))) %counts must have a count for each combo of positive natural symbols and also combinations of overflows (0)
                error('pidgeon hole issue');
            end

        
            obj.total_iterations = uint64(0);
            obj.maxCount = bitcmp(uint64(0));
            obj.nSymbolsPerDim = biggestSymbolPerDimension+1;
            obj.nDims = length(obj.nSymbolsPerDim);

            warning('you need to MANUALLY ensure that prod(1+biggestSymbolPerDimension) is less than 2^64-1 and that no element of counts exceeds 2^64-1')
            
            obj.counts = uint64(counts);
            if(sum(obj.counts<0))
                error('counts must be nonnegative\n');
            end
            obj.counts(obj.counts==0)=uint64(1); 

            obj.perm = 1:(obj.nSymbols);
            [obj.counts,b] = sort(obj.counts,'descend');
            obj.perm = obj.perm(b);
            for idx = 1:obj.nSymbols
                obj.iperm(obj.perm(idx))=idx;
            end
        end
        
        %precondition: current model is valid.
        %e.g. denominator is below near miss.
        function updateModel(obj,symbol,varargin)
            
            obj.total_iterations=obj.total_iterations+1;
            
            symbolIdx = obj.getLinearIdxFromSymbolTuple(symbol);

            if(obj.counts(symbolIdx)>= bitcmp(uint64(0)))
                obj.counts = bitshift(obj.counts-uint64(1),-1);
                obj.counts = obj.counts + uint64(1);
                obj.counts(symbolIdx) = obj.counts(symbolIdx)+uint64(1);
            else
                obj.counts(symbolIdx) = obj.counts(symbolIdx)+1;
            end
            
            [obj.counts,newperm] = sort(obj.counts,'descend');%not efficient but who cares.
            obj.perm = obj.perm(newperm);
            obj.iperm(obj.perm) = 1:obj.nSymbols;

        end
        
        %next step- rewrite sub2ind, ind2sub so that they are more
        %efficient and memoize. 
        function linearIdx = getLinearIdxFromSymbolTuple(obj,symbol)
           dummy = num2cell(symbol+1);
           indexTuple = sub2ind(obj.nSymbolsPerDim, dummy{:});
           linearIdx = obj.iperm(indexTuple);  
        end
        
        function symbol = getSymbolTupleFromLinearIdx(obj,linearIdx)
           linearIdx = obj.perm(linearIdx);
           dummy = cell(1,obj.nDims);
           [dummy{:}] = ind2sub(obj.nSymbolsPerDim,linearIdx);
           indexTuple = cell2mat(dummy);
           symbol = indexTuple-1;
        end
             
             
    end
end 

