classdef simpleKalmanFilterTracking < handle

    properties
        A %system matrix
        W %process noise
        
        C %measurement matrix
        V %measurement noise
        
        xpred %a priori estimator/prediction (t+1|t)
        Ppred% a priori estimator/prediction covariance  (t+1|t)
        
        xpost %posterior estimator/prediction (t|t)
        Ppost %posterior estimator/prediction covariance  (t|t)

        tpred 
        tpost
        
        
    end
    
    properties(Access = private)
        I
    end
    
    methods
        function obj = simpleKalmanFilterTracking(A,C,W,V,xest_init,Pinit_est)
            
            obj.A = A; %system matrix, must be square
            obj.C = C; %measurement matrix
            obj.W = W; %process noise
            obj.V = V; %measurement noise
            
            obj.xpost = xest_init;
            obj.Ppost = Pinit_est;

            obj.xpred = xest_init;
            obj.Ppred = Pinit_est; %trust me this makes sense no matter 
            %if you start by updating or not so long as you call the
            %predict and update methods in the order you want. 
            obj.I = eye(size(A,1));            

            obj.tpred = 1;            
            obj.tpost = 0;
        end
        
        %can only call once per time step
        %not efficient. think about this. 
        %given control input at time t-1, and some measurements from times
        %1....,t, computes best estimate of x_t|(measurements recieved until time t)
        
        function predictUpdate(obj,controlInput)
            obj.xpred = obj.A*obj.xpost;
            obj.Ppred = obj.A*obj.Ppost*obj.A'+obj.W;
            obj.tpred = obj.tpost+1;
        end
                 
        function measurementUpdate(obj,measurement)
            %warning('can only call once per time step');
            innovation = measurement-obj.C*obj.xpred;
            dummy = obj.Ppred*obj.C'; %this memoization might
            %be numerically bad
            S = obj.C*dummy+obj.V;
            Kkf = dummy/S;
            obj.xpost = obj.xpred+Kkf*innovation;
            obj.Ppost = (obj.I-Kkf*obj.C)*obj.Ppred;
            obj.tpost = obj.tpost + 1;
        end
    end
 
end

