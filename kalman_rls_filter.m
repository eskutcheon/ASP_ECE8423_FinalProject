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

function [filter_state, y, e] = kalman_rls_filter(filter_state, x, d)
    % Kalman RLS filter implementation
    TOL = 1e-8;
    M = size(filter_state.F, 1);
    N = size(filter_state.C, 1);
    I_M = eye(M);   % identity matrix for state dimension
    

    % measurement prediction based on current x_hat, of shape (N, 1)
    y_hat = filter_state.C * filter_state.x_hat;
    % innovations vector, of shape (N, 1)
    alpha = d - y_hat;
    % alpha = d - x' * x_hat;
    fprintf("alpha shape: (%d,%d)\n", size(alpha))
    % Innovation covariance - should be (N, N)
    %R = x' * filter_state.K * x + filter_state.Q2 + TOL;
    R = filter_state.C * filter_state.K * filter_state.C' + filter_state.Q2 + TOL * eye(N);
    fprintf("R shape: (%d,%d)\n", size(R))
    % Kalman Gain - should be shape (M, N)
    %G = filter_state.F * filter_state.K * filter_state.C' * pinv(R);
    %G = filter_state.F * filter_state.x_hat * pinv(R);
    G = filter_state.K * filter_state.C' * inv(R);
    fprintf("G shape: (%d,%d)\n", size(G))
    % State Update
    filter_state.x_hat = filter_state.F * filter_state.x_hat + G * alpha;
    fprintf("x_hat shape: (%d,%d)\n", size(filter_state.x_hat))
    % a priori error covariance estimation K(n,n)
    K_priori = filter_state.K - inv(filter_state.F) * G * filter_state.C * filter_state.K;
    %fprintf("K(n,n) shape: (%d,%d)\n", size(G))
    % Covariance Update
    %filter_state.K = filter_state.F * K_priori * filter_state.F' + filter_state.Q1;
    filter_state.K = (I_M - G * filter_state.C) * filter_state.K;
    % Covariance update
    %K_new = filter_state.F * (filter_state.K - G * filter_state.C * filter_state.K) * filter_state.F' + filter_state.Q1;
    fprintf("K(n+1,n) shape: (%d,%d)\n", size(filter_state.K))
    % Compute output and error
    y = filter_state.C * filter_state.x_hat;
    e = d - y;
    %filter_state.K = K_new;
end