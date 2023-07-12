%This function first comptues the integer part of a delta-uniformly 
%quantized version of the argument r. The quantization is an integer, or
%tuple of integers. delta is an option parameter, 1 by default. It then
%threads the integer (or tuple of integers) mapping 0 to 1, 1 to 2, -1 to 3
%, 2  to 4, -2 to 5 and so on, elementwise for each member of the tuple.
%The output, symbol, is a (tuple of) strictly positive natural number(s).

function symbol = quantizeAndThread(r,varargin)
    
    if(isempty(varargin))
        delta=1;
    else
        delta =varargin{1};
    end
    
    symbol = myRound(r./delta);
    symbol = (symbol<=0).*(-2.0*symbol+1)+(symbol>0).*(2.0*symbol);
    
end

function rd = myRound(r)

rd=((r<0).*(mod(2*r,1)==0)).*ceil(r)+(~((r<0).*(mod(2*r,1)==0))).*round(r);

end

