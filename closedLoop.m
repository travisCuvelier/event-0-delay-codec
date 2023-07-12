clear all
close all

addpath('elias_omega')

%inverted pendulum model from https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlDigital

%plant model

Act = [-1.1680 -.0886; 2.003, -.2443]
plantDim = size(Act,1);
Bct =  eye(plantDim)*.2

sigmaInit = .05*eye(plantDim);



Tsample = 1/100;

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
numIterations = 100000;

x = zeros(plantDim,numIterations);
xhat = zeros(plantDim,numIterations);
y = zeros(plantDim,numIterations);
x(:,1) = sigmaInit*randn(plantDim,1);
xpred_init = zeros(plantDim,1); %initial plant estimate
Ppred_init = sigmaInit^2*eye(plantDim);%"initial" post measurement KF 



%set up encoder, decoder, quantizer, etc

cutoffs = [127,127]; %if cutoff= [21,21,21,21] the quantizer will 
%have bins from -10 to 10 in intevals of 1 unit. Anything outside that must
%be sent another way (elias)
%rate distortion estimate is 13.5 bits or about 3.3 
%bits per dimension. 2^3.3 = 10ish. We doubled it.
    
uniformPrior = uint64(ones(prod(cutoffs+1),1));
delta = eye(size(Act,1));%sqrtm(12}{":"*policyStruct.V); to save time I hard coded delta
herald = 0; %don't change this


dither = rand(plantDim,numIterations)-.5;

encoderModel = sortedAdaptiveCutoffPMF64(uniformPrior,cutoffs);
encoder = eventEncoder64(encoderModel);

Dct = 4e-4; %continuous-time distortion contraint
Ddt = Tsample*Dct-bline;
if(Ddt<0)
    error('sampling rate needs to be higher (Tsample smaller) to sustain this performance')
end

sensingPolicy = rateDistortionTracking(A,W,Q,Ddt,'mosek'); %about 13.5 bits
C = sensingPolicy.C; 
V = sensingPolicy.V;

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
    symbols = quantizeAndThread(measurement+dither(:,iteration)).'; %asssumes delta = eye(4)

    %figure our the measurement the decoder will eventually recieve (either
    %this iteration or the next iteration) and update the encoder's kalman
    %filter accordingly
    measurementDecoderEventuallyReceives = unthreadAndReconstruct(symbols).'+C*encoderKF.xpred-dither(:,iteration);
    encoderKF.measurementUpdate(measurementDecoderEventuallyReceives);
    encoderKF.predictUpdate();

    %figure out if there are any overflows, and if so save them to encode
    %next time. 
    txOverflows = symbols(symbols>cutoffs);
    txNumOverflows = numel(txOverflows);
    totalOverflows = totalOverflows+txNumOverflows; %count overflows over simulation


    %replace "overflows" with escape symbol (0). 
    symbolsToEncodeNormally = symbols;
    symbolsToEncodeNormally(symbolsToEncodeNormally>cutoffs) = herald; %escape symbols
    %the herald indicates that the beginning of the next tramsission will
    %consist of overflowed symbols
    
    %encode the tuple with overflows replaced by escape using the
    %nonsingular code. update probability model at encoder accordingly
    normalBits = encoder.encodeSymbol(symbolsToEncodeNormally);
    encoderModel.updateModel(symbolsToEncodeNormally) %models only depend on "normally encoded bits" so they are updated each time regardless.

    %keep track of the bits count not including the overflows
    cwlNoOverflows(iteration)=length(normalBits);

    bits = [bits, normalBits];
    cwl(iteration) = numel(bits);

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

    tcost(iteration) = (x(:,iteration)-estimates(:,iteration))'*Q*(x(:,iteration)- estimates(:,iteration));
    x(:,iteration+1) = A*x(:,iteration)+sqrtm(W)*randn(plantDim,1);


end

