function out = omegaDecodeStreaming(binArray)

    out = 1; 
    N=-1;
    while(binArray(N+2)~=0)
        nplus1word = binArray(N+2:N+2+out);
        N = N+out+1; 
        out = bi2de(nplus1word,'left-msb'); 
    end
  

end