Act = [-1.1680 -.0886; 2.003, -.2443] %put your A and B here. 
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

Dct = .01; %% THIS IS WHAT YOU CHANGE; your continuous time distortion. 

Ddt = Tsample*Dct-bline;

sensingPolicy = rateDistortionTracking(A,W,Q,Ddt,'mosek'); 
dtLowerBound = theta_inv(sensingPolicy.minimumBits)/Tsample