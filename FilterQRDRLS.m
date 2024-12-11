classdef FilterQRDRLS < RLSVariantBase
    properties (Access = protected)
        P; % Square-root covariance matrix
    end
    
    methods
        function initialize(obj, filterOrder, lambda, delta)
            obj.lambda = lambda;
            obj.delta = delta;
            obj.initializeWeights(filterOrder);
            obj.P = delta * eye(filterOrder);
        end
        
        function updateMatrix(obj, input)
            % Update square-root covariance matrix using QR decomposition
            obj.P = obj.lambda^(-0.5) * qr([obj.P, input], 0);
        end
        
        function gain = computeGain(obj, input)
            % Compute gain vector using square-root covariance matrix
            gain = obj.P \ (input / (obj.lambda + input' * (obj.P' * obj.P) * input));
        end
    end
end
