function bits = omegaEncode(integer)
    
    if(integer<=0 || integer~=floor(integer))
        error('The omega encoder only accepts positive integer inputs');
    end
    
    bits = [dive(integer),0];
    
end

function bits = dive(N)

    if(N ==1)
        bits = [];
    else
        st = de2bi(N,'left-msb'); %big endian is important
        bits = [dive(numel(st)-1),st];
    end

end