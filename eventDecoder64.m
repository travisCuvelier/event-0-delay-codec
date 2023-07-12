classdef eventDecoder64 < handle
   
    properties(Constant)
        precision = 64;
    end
    
    properties (SetAccess = public, GetAccess = public)
        model
        nSymbols
        %This script does not update the model. 
    end

    properties (SetAccess = private, GetAccess = public)
        codebook
        %This script does not update the model. 
    end
    
    methods(Static=true)
        
        %big endian binary expansion of z (assumed column) given bytewidth.
        function ze = binaryExpand(z,bytewidth)
            ze = zeros(numel(z),bytewidth);
            
            for idx = 1:bytewidth
                ze(:,bytewidth+1-idx) = mod(z,2);
                z = z - ze(:,bytewidth+1-idx);
                z = bitshift(z,-1);
            end
        end
        
        %this is like bi2de(,'left-msb') in the communications toolbox.
        function d = binary2decimal(b)
            
            [~,nc] = size(b);
            
            d = sum(uint64(2.^((nc-1):-1:0)).*uint64(b),2);
            
        end %the rate of the code in bits/channel use
        
    end
    
    methods
        %constructor
        function obj = eventDecoder64(model)
            obj.model  = model;   
            obj.nSymbols = model.nSymbols;
            obj.codebook = cell(obj.nSymbols,1);
            obj.codebook{1} = [];
            nindicesprev = uint64(1);
            nbits = uint64(1);
            nindiceshere = uint64(2^nbits)+uint64(nindicesprev);
         
            for idx = 2:obj.nSymbols
               if(idx > nindiceshere)
                   nindicesprev = nindiceshere;
                   nbits = uint64(nbits)+uint64(1);
                   nindiceshere =nindicesprev+uint64(2^nbits);
               end
                z = uint64(idx -1-nindicesprev);
                obj.codebook{idx} = eventEncoder64.binaryExpand(z,nbits);
            end
            
        end
        
        function symbol = decodeSymbol(obj,codeword)

            nbits = numel(codeword);
            if(nbits>obj.precision)
                error('this shouldn''t happen. codeword length has to be (much) less than precision')
            end

            if(numel(codeword)==0)
                linearIndex = 1;
            else

                base = uint64(1);
                for idx = 1:(nbits-1)
                    base = uint64(base) + uint64(2^idx);
                end 
    
                linearIndex = base+1+eventDecoder64.binary2decimal(codeword);
    
            end

            symbol = obj.model.getSymbolTupleFromLinearIdx(linearIndex);

        end        
            
    end
    
    
end