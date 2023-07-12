% add the following directories to your MATLAB path
%see notes in YALMIP-master
% ->/YALMIP-master
% ->/YALMIP-master/extras
% ->/YALMIP-master/solvers
% ->/YALMIP-master/modules
% ->/YALMIP-master/modules/parametric
% ->/YALMIP-master/modules/moment
% ->/YALMIP-master/modules/global
% ->/YALMIP-master/modules/sos
% ->/YALMIP-master/operators
% addpath('../../YALMIP-master')
% addpath('../../YALMIP-master/extras')
% addpath('../../YALMIP-master/solvers')
% addpath('../../YALMIP-master/modules')
% addpath('../../YALMIP-master/modules/parametric')
% addpath('../../YALMIP-master/modules/moment')
% addpath('../../YALMIP-master/modules/global')
% addpath('../../YALMIP-master/modules/sos')
% addpath('../../YALMIP-master/operators')

% addpath('YALMIP-master')
% addpath('YALMIP-master/extras')
% addpath('YALMIP-master/solvers')
% addpath('YALMIP-master/modules')
% addpath('YALMIP-master/modules/parametric')
% addpath('YALMIP-master/modules/moment')
% addpath('YALMIP-master/modules/global')
% addpath('YALMIP-master/modules/sos')
% addpath('YALMIP-master/operators')

%varargin is a string that specifies solver, e.g. 'mosek' for MOSEK,
%'sdpt3' for SDPT3...

%Make sure you have yalmip ready
%to go and solvers appropriately working before you call this. 

%Asys is the SAMPLED system feedback matrix, Wcov is the SAMPLED
%process noise COVARIANCE (BB' in our paper), 
%Qsemidef is a PSD matrix, same dimensions as A, that weighs distortion. 
%Dpos is a positive distrotion target.
%reccomend calling with optional argument 'mosek' if you have it installed
%this function doesn't check for most input errors. Make sure A and W are
%sampled at the same rate, distrotion, Q sampling rate appropriate. 
%Make sure your distortion constraint is feasible, not vacuous, etc. 


function policy = rateDistortionTracking(Asys,Wcov, Qsemidef,Dpos,varargin)

    if(isempty(varargin))
        solver = 'sdpt3';
    else
        solver = varargin{1};
    end
    

    if(rank(Wcov)~=size(Wcov,1))
        warning('your process noise covariance is supposed to be full rank. proceed at your own peril.')
    end
    
    n = size(Asys,1);
    P = sdpvar(n);
    PI = sdpvar(n);
    F = [ P>=0, PI>=0,trace(Qsemidef*P)<= Dpos,Asys*P*Asys'+Wcov-P>=0, [P-PI,P*Asys';Asys*P,Asys*P*Asys'+Wcov]>=0  ];
    
    diagnostic=optimize(F,-logdet(PI),sdpsettings('solver',solver,'verbose',0));
    if(diagnostic.problem)
        s = sprintf('%d\n',diagnostic.problem);
        warning(strcat('issues with solver, error code = ',s))
        diagnostic
    end

    
    policy.minbitsa = .5*( log2(det(Asys*value(P)*Asys'+Wcov))- log2(det(value(P))) );
    policy.minbitsb = .5*(-log2(det(value(PI))) +log2(det(Wcov)));
    policy.minimumBits = .5*(-log2(det(value(PI))) +log2(det(Wcov)));
    policy.P = value(P);
    policy.PI = value(PI);
    policy.SNR = inv(policy.P)-inv(Asys*policy.P*Asys'+Wcov);
    [U,E,~] = svd(policy.SNR);
    rk = 0;
    rplus = 1;
    while(rk < size(E,1) && E(rplus,rplus) ~=0)
        rk = rplus;
        rplus = rplus+1;
    end
    policy.rank = rk;
    Uprime = U*sqrt(E)/sqrt(12);
    policy.C = (Uprime(:,1:rk))'; %sensor gain
    policy.V = eye(rk)/12; %sensor noise, set to use uniform dithered
    % quantization on integer cells. 
    policy.solverDiagnostics = diagnostic;
    
end