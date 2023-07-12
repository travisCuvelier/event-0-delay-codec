function recon = unthreadAndReconstruct(n,varargin)
    
    if(isempty(varargin))
        delta=1;
    else
        delta =varargin{1};
    end
    
    recon = (mod(n,2)==0).*(delta.*n)/2-(mod(n,2)).*(delta.*(n-1))/2;
        
end

