classdef LandmarkRangeBearingEdge < g2o.core.BaseBinaryEdge
    % LandmarkRangeBearingEdge summary of LandmarkRangeBearingEdge
    %
    % This class stores an edge which represents the factor for observing
    % the range and bearing of a landmark from the vehicle. Note that the
    % sensor is fixed to the platform.
    %
    % The measurement model is
    %
    %    z_(k+1)=h[x_(k+1)]+w_(k+1)
    %
    % The measurements are r_(k+1) and beta_(k+1) and are given as follows.
    % The sensor is at (lx, ly).
    %
    %    dx = lx - x_(k+1); dy = ly - y_(k+1)
    %
    %    r(k+1) = sqrt(dx^2+dy^2)
    %    beta(k+1) = atan2(dy, dx) - theta_(k+1)
    %
    % The error term
    %    e(x,z) = z(k+1) - h[x(k+1)]
    %
    % However, remember that angle wrapping is required, so you will need
    % to handle this appropriately in compute error.
    %
    % Note this requires estimates from two vertices - x_(k+1) and l_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k+1) and slot 2 contains l_(k+1).
    
    methods(Access = public)
    
        function obj = LandmarkRangeBearingEdge()
            % LandmarkRangeBearingEdge for LandmarkRangeBearingEdge
            %
            % Syntax:
            %   obj = LandmarkRangeBearingEdge();
            %
            % Description:
            %   Creates an instance of the LandmarkRangeBearingEdge object.
            %   Note we feed in to the constructor the landmark position.
            %   This is to show there is another way to implement this
            %   functionality from the range bearing edge from activity 3.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a LandmarkRangeBearingEdge.

            obj = obj@g2o.core.BaseBinaryEdge(2);
        end
        
        function initialEstimate(obj)
            % INITIALESTIMATE Compute the initial estimate of the landmark.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial estimate of the landmark given the
            %   platform pose and observation.

            x = obj.edgeVertices{1}.x;
            r = obj.z(1);
            beta = obj.z(2);

            phi = x(3) + beta;
            lxy = x(1:2) + r * [cos(phi); sin(phi)];
            obj.edgeVertices{2}.setEstimate(lxy);
        end
        
        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the predicted and actual range-bearing measurement.

            x = obj.edgeVertices{1}.estimate();
            l = obj.edgeVertices{2}.estimate();

            dx = l(1) - x(1);
            dy = l(2) - x(2);
            r = sqrt(dx^2 + dy^2);
            h = [r; atan2(dy, dx) - x(3)];

            obj.errorZ = obj.z - h;
            obj.errorZ(2) = g2o.stuff.normalize_theta(obj.errorZ(2));
        end
        
        function linearizeOplus(obj)
            % linearizeOplus Compute the Jacobian of the error in the edge.
            %
            % Syntax:
            %   obj.linearizeOplus();
            %
            % Description:
            %   Compute the Jacobian of the error function with respect to
            %   the vertex.
            %

            x = obj.edgeVertices{1}.x;
            l = obj.edgeVertices{2}.x;

            dx = l(1) - x(1);
            dy = l(2) - x(2);
            r2 = dx^2 + dy^2;
            r = sqrt(r2);

            % Jacobian wrt platform state x=[x;y;theta]
            obj.J{1} = [dx/r dy/r 0; ...
                -dy/r2 dx/r2 1];

            % Jacobian wrt landmark state l=[lx;ly]
            obj.J{2} = [-dx/r -dy/r; ...
                dy/r2 -dx/r2];
        end        
    end
end
