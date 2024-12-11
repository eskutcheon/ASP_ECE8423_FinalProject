classdef FilterBasicRLS < RLSFilterBase
    properties (Access = protected)
        w;          % Filter weights
        eps;        % weight update error
        %err_init;  % initial error value to use in one formulation of the error that I tested
        P;          % Inverse covariance matrix
        lambda;     % Forgetting factor
        delta;      % Regularization parameter
    end

    methods
        function obj = initialize(obj, filter_order, lambda, delta)
            % Initialize filter parameters
            obj.w = zeros(filter_order, 1);             % initialize weights to zero
            obj.eps = ones(filter_order, 1);            % weight update error initialized to ones
            %obj.err_init = obj.eps(1);
            obj.P = (1 / delta) * eye(filter_order);    % regularized identity matrix
            obj.lambda = lambda;
            obj.delta = delta;
        end

        function [obj, y, err] = adapt(obj, x, d)
            % gain vector (Kalman gain)
            k = obj.compute_gain(x);
            % a priori estimation error
            xi = obj.compute_priori_error(x, d);
            % a posteriori weight error (previous error - k*xi)
            %obj.eps = obj.eps - k * xi';
            % update weights with a priori error
            obj.w = obj.w + k * xi';
            % output signal (filtered signal)
            y = obj.compute_output(x);
            % a posteriori estimation error (desired - output)
            err = obj.compute_posteriori_error(y, d);
            %err = obj.err_init - sum(k .* xi);
            % seemingly equivalent definition with the conversion factor
            %err = (1 - k' * x) * xi;
            %gamma = err / xi;
            % update inverse covariance matrix
            obj.P = (1 / obj.lambda) * (obj.P - k * x' * obj.P);
        end

        % --- Helper Functions with Conventional Variables ---
        function k = compute_gain(obj, x)
            TOL = 1e-8;
            %%
            % Compute Kalman gain - simplified form without $\lambda^{-1}$
            k = (obj.P * x) / (obj.lambda + x' * obj.P * x + TOL);
        end
    end
end
