clear all
close all

addpath('elias_omega')

%plant model
%https://www.informit.com/articles/article.aspx?p=32090&seqNum=2 operating
%condition one from Exothermic CSTR


Act = [-1.1680 -.0886; 2.003, -.2443]
plantDim = size(Act,1);
Bct =  eye(plantDim)*.2;

sigmaInit = .05*eye(plantDim);

Tsample = 1/100;

%Build discrete time model
syms x
f(x) =  expm(Act*x)*(Bct*Bct')*expm(Act'*x);
g = matlabFunction(int(f));
fun2 = @(t) g(t)-g(0);
bline = trace(integral(fun2,0,Tsample,'ArrayValued',true));
clear x
clear f

A = expm(Act*Tsample);
W = integral(@(x) expm(Act*x)*(Bct*Bct')*expm(Act'*x),0,Tsample,'ArrayValued',true);
Q = integral(@(x) expm(Act*x)*expm(Act'*x),0,Tsample,'ArrayValued',true);




%set up some simulation variables
numIterations = 10000;


%if you want to loop over distortions, start here. can move somethings out
%of the loop like "cuttoffs" and "uniformPrior" if you're not going to
%change the cuttoffs at each distortion. I didn't do this due the thematic
%organization

%more simulation setup 
x = zeros(plantDim,numIterations);
xhat = zeros(plantDim,numIterations);
y = zeros(plantDim,numIterations);
x(:,1) = sigmaInit*randn(plantDim,1);
xpred_init = zeros(plantDim,1); %initial plant estimate
Ppred_init = sigmaInit^2*eye(plantDim);%"initial" post measurement KF 



%set up encoder, decoder, quantizer, etc

cutoffs = [127,127]; %if cutoff= [21,21], the quantizer will 
%have bins from -10 to 10 in intevals of 1 unit. Dont make
%(prod(cuttoffs+1) anywhere near 2^64-1).
    
uniformPrior = uint64(ones(prod(cutoffs+1),1));

herald = 0; %don't change this
dither = rand(plantDim,numIterations)-.5;%dither for quantization to 
                                         %integer lattice. don't change

encoderModel = sortedAdaptiveCutoffPMF64(uniformPrior,cutoffs);
encoder = eventEncoder64(encoderModel);

Dct = 4e-4; %continuous-time distortion contraint. this is what you'd 
            %want to change in your loop. you might also want to change
            %cuttoffs and/or sampling rate. 

Ddt = Tsample*Dct-bline; %discrete-time distortion constraint
if(Ddt<0)
    error('sampling rate needs to be higher (Tsample smaller) to sustain this performance')
end

%solve the rate-distortion problem
sensingPolicy = rateDistortionTracking(A,W,Q,Ddt,'mosek'); 
C = sensingPolicy.C; %fig 2 in the paper explains what these are.
V = sensingPolicy.V;
%sensingPolicy.minimumBits is R(D_{d,\tau}, \overline{Q}_{A.\tau},\tau) 
%from the paper 

%discrete time lower bound calculation
dtLowerBound = theta_inv(sensingPolicy.minimumBits)/Tsample;

encoderKF = simpleKalmanFilterTracking(A,C,W,V,xpred_init,Ppred_init);

decoderModel = sortedAdaptiveCutoffPMF64(uniformPrior,cutoffs);
decoder = eventDecoder64(decoderModel);
decoderKF = simpleKalmanFilterTracking(A,C,W,V,xpred_init,Ppred_init);
totalOverflows = 0;
rxNumOverflows = 0;
txNumOverflows = 0;
estimates = zeros(plantDim,numIterations);

for iteration = 1:numIterations



    if(mod(iteration,1000)==0)
        fprintf('%3d %% iteration of CT Cost %5.2f \n',iteration*100/numIterations,Dct)
    end


    %%The encoder

    bits = [];

    %first, encode overflow bits from the last iteration
    if(txNumOverflows ~= 0) 
        for overIdx = 1:txNumOverflows
            bits = [bits, omegaEncode(txOverflows(overIdx))];
        end
    end
    
    %make measurement with respect to the encoder's Kalman filter.
    measurement = C*(x(:,iteration)-encoderKF.xpred);

    %turn measurement into a tuple of positive integers. 
    symbols = quantizeAndThread(measurement+dither(:,iteration)).'; 
    %quantizeAndThread first quantizes measurement to integer lattice, 
    %then bijectivelly wraps these points to a tuple of positive integers.
    %this tuple is easier to encode with this cutoff method. pretty much,
    %symbols(idx) overflows if symbols(idx)>cutoffs(idx). 


    %the encoder figures out the measurement the decoder will eventually
    %recieve (either this iteration or the next iteration depending on
    %whether or not there was an overflow) and update its KF accordingly. 
    measurementDecoderEventuallyReceives = ...
        unthreadAndReconstruct(symbols).'+C*encoderKF.xpred-...
                                                    dither(:,iteration);

    encoderKF.measurementUpdate(measurementDecoderEventuallyReceives);
    encoderKF.predictUpdate();


    %figure out if there are any overflows, and if so save them to encode
    %next time. 
    txOverflows = symbols(symbols>cutoffs);
    txNumOverflows = numel(txOverflows);
    totalOverflows = totalOverflows+txNumOverflows; %track total overflows
                                                    %over simulation. if
                                                    %you loop this over
                                                    %multiple distortions,
                                                    %maybe index this
                                                    %variable to keep track
                                                    %of overflows at each
                                                    %distortion level. 

    %replace "overflows" with escape symbol (0). 
    symbolsToEncodeNormally = symbols;
    symbolsToEncodeNormally(symbolsToEncodeNormally>cutoffs) = herald; 
    %herald = the escape symbol
    %the herald indicates that the beginning of the next tramsission will
    %consist of overflowed symbols
    
    %encode the tuple with overflows replaced by escapes using the
    %nonsingular code. update probability model at encoder accordingly
    normalBits = encoder.encodeSymbol(symbolsToEncodeNormally);
    encoderModel.updateModel(symbolsToEncodeNormally) %models only depend 
                                                      %on "normally
                                                      %encoded bits" so
                                                      %they are updated 
                                                      %every time.


    %keep track of the bits count not including the overflows
    cwlNoOverflows(iteration)=length(normalBits);

    bits = [bits, normalBits]; %total bits is overflow encoding from 
                               %beginning plus this time's NS encoding. 

    cwl(iteration) = numel(bits); %track total num bits. If you loop this 
                                  %over distrotions, you'll want to 
                                  %compute the mean of this for each
                                  %distortion level. 
                                  
    %% the decoder

    %if on the previous iteration the recieved transmission included
    %escape symbols (encoded with the nonsingular code), then this
    %iteration the encoder will send the real value of those overflowed
    %symbols using the Elias code at the beginning of its transmission. The
    %deocder knows how many overflows to expect because it counted the
    %number of heralds recieved last time. 
    if(rxNumOverflows~=0)  
  
        %decode each overflowed symbol. Since the Elias code is prefix, we
        %remove the bits used to encode the overflow symbols sequentially
        %(ie, we know when each overflow symbol's codeword ends via the
        %prefix property)
        rxOverflows = zeros(rxNumOverflows,1);
        for rxOverIdx = 1:rxNumOverflows
            rxOverflows(rxOverIdx) = omegaDecode(bits);
            numBitsToDrop = numel(omegaEncode(rxOverflows(rxOverIdx)));
            bits = bits(numBitsToDrop+1:end);
        end

        %these are still the symbols from last time. Replace the heralds
        %with their true overflow values. 
        rxSymbols(rxSymbols == herald) = rxOverflows;

        %process measurement, update Kalman filter.  
        rxMeasurement = unthreadAndReconstruct(rxSymbols).'+C*decoderKF.xpred-dither(:,iteration-1);%need to use LAST TIME's dither, since we are fixing last time's tranmsission
        decoderKF.measurementUpdate(rxMeasurement);
        decoderKF.predictUpdate(); %the encoder and decoder KFs are now 
                                   %synched up again.         
        
        rxNumOverflows = 0;
        rxOverflows = [];

        %we have now fixed our previous transmission. Now, the encoder and
        %decoder KF's should be synchronized. remaining "bits" were encoded
        %with the nonsingular code. They do not contain encoded overflows
        %but may include "herald" symbols.
    end


    %decode the part of the codeword that was encoded with the nonsingular
    %code, and update the decoder's model as such. The encoder and decoder
    %models are always synchronized, regardless of whether or not there
    %were overflows. 
    rxSymbols = decoder.decodeSymbol(bits);
    decoderModel.updateModel(rxSymbols);

    %if heralds were found in *this* nonsingular codeword, next time the
    %decoder needs to expect a prefix containing the overflow symbols
    %encoded with the Elias universal code. It counts the herald to
    %recognize how many to expect.

    rxNumOverflows = sum(rxSymbols==herald);
    
    if(rxNumOverflows ~= 0)
        %no measurments this time, so the prediction update is the best we
        %have. Encoder and decoder KF will be out of synch until the
        %overflow symbols are recovered. 

        estimates(:,iteration) = decoderKF.xpred; %we didn't recieve a measurement this time.

    else %if no overflows, everything is great
        %compute the measurement and associated kalman filter estimate.  
        rxMeasurement = unthreadAndReconstruct(rxSymbols).'+C*decoderKF.xpred-dither(:,iteration);%use this time's dither

       
        decoderKF.measurementUpdate(rxMeasurement); %KF measurement update
        estimates(:,iteration) = decoderKF.xpost; %posterior estimate is the best we have
        decoderKF.predictUpdate(); %run prediction. 

        %if there are overflows to be recieved in the next iteration, this
        %block will not be run. The measurement and predict updates for the
        %decoder KF will be run at the beginning of the next iteration in
        %the if(rxNumOverflows~=0)  block.
    end

    %compute tracking error. Compare this to Ddt. 
    tcost(iteration) = (x(:,iteration)-estimates(:,iteration))'*Q*...
                                (x(:,iteration)- estimates(:,iteration));

    %update plant. 
    x(:,iteration+1) = A*x(:,iteration)+sqrtm(W)*randn(plantDim,1);


end

commCostBitsPerSec = mean(cwl)/Tsample; %bits per second comm cost. 

