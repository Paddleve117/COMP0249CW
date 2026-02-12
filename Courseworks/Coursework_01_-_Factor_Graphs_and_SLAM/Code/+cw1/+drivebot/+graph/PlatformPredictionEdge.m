classdef PlatformPredictionEdge < g2o.core.BaseBinaryEdge
    % PlatformPredictionEdge summary of PlatformPredictionEdge
    %
    % This class stores the factor representing the process model which
    % transforms the state from timestep k to k+1
    %
    % The process model is as follows.
    %
    % Define the rotation vector
    %
    %   M = dT * [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0;0 0 1];
    %
    % The new state is predicted from 
    %
    %   x_(k+1) = x_(k) + M * [vx;vy;theta]
    %
    % Note in this case the measurement is actually the mean of the process
    % noise. It has a value of 0. The error vector is given by
    %
    % e(x,z) = inv(M) * (x_(k+1) - x_(k))
    %
    % Note this requires estimates from two vertices - x_(k) and x_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k) and slot 2 contains x_(k+1).
    
    properties(Access = protected)
        % The length of the time step
        dT;
    end
    
    methods(Access = public)
        function obj = PlatformPredictionEdge(dT)
            % PlatformPredictionEdge for PlatformPredictionEdge
            %
            % Syntax:
            %   obj = PlatformPredictionEdge(dT);
            %
            % Description:
            %   Creates an instance of the PlatformPredictionEdge object.
            %   This predicts the state from one timestep to the next. The
            %   length of the prediction interval is dT.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a PlatformPredictionEdge

            assert(dT >= 0);
            obj = obj@g2o.core.BaseBinaryEdge(3);            
            obj.dT = dT;
        end
       
        function initialEstimate(obj)
            % INITIALESTIMATE Compute the initial estimate of a platform.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial estimate of the platform x_(k+1) given
            %   an estimate of the platform at time x_(k) and the control
            %   input u_(k+1)

            xk = obj.edgeVertices{1}.x;
            theta = xk(3);
            M = obj.dT * [cos(theta) -sin(theta) 0; ...
                sin(theta) cos(theta) 0; ...
                0 0 1];

            % Compute the posterior assuming no process noise
            xkp1 = xk + M * obj.z;
            xkp1(3) = g2o.stuff.normalize_theta(xkp1(3));
            obj.edgeVertices{2}.x = xkp1;
        end
        
        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the measurement and the parameter state in the
            %   vertex. Note the error enters in a nonlinear manner, so the
            %   equation has to be rearranged to make the error the subject
            %   of the formulat
                       
            xk = obj.edgeVertices{1}.x;
            xkp1 = obj.edgeVertices{2}.x;

            theta = xk(3);
            M = obj.dT * [cos(theta) -sin(theta) 0; ...
                sin(theta) cos(theta) 0; ...
                0 0 1];

            obj.errorZ = M \ (xkp1 - xk) - obj.z;
            obj.errorZ(3) = g2o.stuff.normalize_theta(obj.errorZ(3));
        end
        
        % Compute the Jacobians
        function linearizeOplus(obj)
            % LINEARIZEOPLUS Compute the Jacobians for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the Jacobians for the edge. Since we have two
            %   vertices which contribute to the edge, the Jacobians with
            %   respect to both of them must be computed.
            %

            xk = obj.edgeVertices{1}.x;
            xkp1 = obj.edgeVertices{2}.x;

            theta = xk(3);
            invM = (1 / obj.dT) * [cos(theta) sin(theta) 0; ...
                -sin(theta) cos(theta) 0; ...
                0 0 1];

            deltaX = xkp1 - xk;
            dinvMdtheta = (1 / obj.dT) * [-sin(theta) cos(theta) 0; ...
                -cos(theta) -sin(theta) 0; ...
                0 0 0];

            obj.J{1} = -invM;
            obj.J{1}(:, 3) = obj.J{1}(:, 3) + dinvMdtheta * deltaX;

            obj.J{2} = invM;
        end
    end    
end
