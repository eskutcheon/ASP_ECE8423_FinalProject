classdef (Abstract) RLSFilterBase
    properties (Abstract, Access = protected)
        w;          % Filter weights
        %eps;        % weight update error
        %P;          % Inverse covariance matrix
        lambda;     % Forgetting factor
        %delta;      % Regularization parameter
    end

    methods (Abstract)
        % abstract initialization function
        %obj = initialize(obj, filter_order, lambda, delta);
        obj = initialize(obj, filter_order, varargin);
        % abstract adaptation function
        [obj, y, err] = adapt(obj, x, d);
    end
    methods
        % --- Helper Functions with Conventional Variables ---
        function y = compute_output(obj, x)
            % Compute filter output
            y = obj.w' * x;
        end
        function xi = compute_priori_error(obj, x, d)
            xi = d - obj.w' * x;
        end
        function err = compute_posteriori_error(obj, y, d)
            % Compute error signal
            err = d - y;
        end
        % --- end helper functions ---
    end
end
