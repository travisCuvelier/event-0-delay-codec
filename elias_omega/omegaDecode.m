function N = omegaDecode(binArray)
    N = omegaDive(binArray,1);
end

function N= omegaDive(binArray,N)

    if(binArray(1) ~=0)
        nplus1word = binArray(1:(N+1));
        N = omegaDive(binArray((N+2):end),bi2de(nplus1word,'left-msb')); 
        %our use of end is suspicious.
    end

end