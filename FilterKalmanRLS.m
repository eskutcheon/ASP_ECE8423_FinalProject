classdef FilterKalmanRLS < RLSFilterBase
    properties (Access = protected)
        w;          % Filter weights
        eps;        % weight update error
        P;          % Inverse covariance matrix
        lambda;     % Forgetting factor
    end
    properties
        F;       % State transition matrix
        C;       % Observation matrix
        Q1;      % Process noise covariance
        Q2;      % measurement noise covariance
    end

    methods
        function obj = initialize(obj, filter_order, varargin)
            % parsing optional arguments
            p = inputParser;
            addOptional(p, 'K', eye(filter_order)); % should be initialized by the NoiseCanceller
            addOptional(p, 'F', eye(filter_order));
            addOptional(p, 'Q1', 1e-5 * eye(filter_order));
            addOptional(p, 'Q2', 1e-2 * eye(filter_order));
            parse(p, varargin{:});
            % Initialize filter parameters
            obj.w = zeros(filter_order, 1);     % initialize weights to zero
            obj.eps = ones(filter_order, 1);    % weight update error initialized to ones
            obj.P = eye(filter_order);          % initial error covariance matrix
            obj.lambda = 1;                     % No forgetting in standard Kalman filter
            % Kalman Filter initializations
            obj.C = eye(filter_order);  % Observation matrix
            obj.F = p.Results.F;
            obj.Q1 = p.Results.Q1;
            obj.Q2 = p.Results.Q2;
        end

        function [obj, y, err] = adapt(obj, x, d)
            % prediction step
            w_hat = obj.F * obj.w;                 % State prediction
            P_hat = obj.F * obj.P * obj.F' + obj.Q1; % Covariance prediction
            fprintf("P_hat shape: (%d,%d)\n", size(P_hat))
            % measurement matrix
            C = x; % Since the observation is scalar, C is the input vector
            fprintf("measurement matrix shape: (%d,%d)\n", size(C))
            % Innovation
            d_hat = C' * w_hat;
            innovation = d - d_hat;
            % Innovation covariance
            R = C' * P_hat * C + obj.Q2;
            fprintf("innovation covariance shape: (%d,%d)\n", size(R))
            % Kalman Gain
            G = (P_hat * C) * inv(R);
            % State Update
            obj.w = w_hat + G * innovation;
            % Covariance Update
            obj.P = (eye(length(obj.w)) - G * C') * P_hat;
            % Return filter output and estimation error
            y = obj.w' * x;
            err = d - y;
        end
    end
end
