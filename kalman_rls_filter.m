%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman filter implementation for the adaptive noise canceller, following
% a state-space model which extends the recursive least-squares (RLS)
% filter with system (state update) equation
%       $x(n+1) = F(n+1,n)x(n) + v_1(n)$,
% where $x(n)$ is the filter state at time $n$, $F(n+1,n)$ is the state
% transition matrix from $n$ to $n+1$, and $v_1$ is zero-mean white noise,
% and measurement (output) equation
%       $y(n) = C(n)x(n) + v_2(n)$,
% where $y(n)$ is the filter output, $C(n)$ is the measurement matrix at
% time $n$, and $v_2$ is zero-mean white noise independent of $v_1$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filter_state, y, e] = kalman_rls_filter(action, M, N, lambda, delta, varargin)
    % Kalman RLS filter implementation
    persistent x K Q1 Q2 F C;
    % NOTE: M = filter length; will use N for the number of samples from here on
    switch action
        case 'initialize'
            % Initialize filter parameters
            if length(varargin) == 0
                x = zeros(M, 1);  % filter state x(n)
                K = eye(M);       % correlation matrix of state prediction error
            else
                x = varargin{1};
                % correlation matrix of state prediction error, initialized with initial state
                x_white = x - mean(x);
                % NOTE: error covariance should be of shape (M, M)
                K = mean(x_white * x_white');
            end
            % process noise variance to form Q1
            v1 = varargin{};
            Q1 = 1e-5 * eye(M);  % covariance matrix of v_1 in the system equation
            Q2 = 1e-2 * eye(N);  % covariance matrix of v_2 in the measurement equation
            F = (1 - delta) .* eye(M);          % state transition matrix
            % FIXME: needs to be N x M
            %C = eye(M);          % FIXME: measurement equation matrix
            C = sprand(N, M, 0.01)          % FIXME: measurement equation matrix
            fprintf("testing the randomly-generated C: %d\n", nnz(C));
            filter_state = struct('x', x, 'K', K, 'Q1', Q1, 'Q2', Q2, 'F', F, 'C', C);

        case 'adapt'
            % initialization step from the noise canceller
            filter_state = varargin{1};
            x = varargin{2};
            d = varargin{3};
            % a priori prediction step
            x_hat = filter_state.F * filter_state.x;
            P_hat = filter_state.F * filter_state.K * filter_state.F' + filter_state.Q1;
            % Innovation - should be shape (N, 1)
            innovation = d - x' * x_hat;
            % Innovation covariance - should be (N, N)
            R = x' * P_hat * x + filter_state.Q2;
            % Kalman Gain - should be shape (M, N)
            G = (P_hat * x) ./ R;

            % State Update
            filter_state.x = x_hat + G * innovation;
            % Covariance Update
            filter_state.K = P_hat - G * x' * P_hat;

            % Compute output and error
            y = filter_state.x' * x;
            e = d - y;
            filter_state = struct('x', filter_state.x, 'K', filter_state.K, ...
                'Q1', filter_state.Q1, 'Q2', filter_state.Q2, 'F', ...
                filter_state.F, 'C', filter_state.C);
    end
end