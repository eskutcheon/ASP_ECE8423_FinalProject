function [filter_state, y, e] = inv_qrd_rls_filter(filter_state, u, d)
    % Inverse QRD-RLS algorithm
    Phi_inv = filter_state.Phi_inv;   % Inverse upper triangular R matrix
    w = filter_state.w;             % Filter weights
    lambda = filter_state.lambda;

    % Compute the Kalman gain vector
    g = Phi_inv * u; % Shape: (M x 1)

    % Compute filter output
    y = w' * u; % Shape: scalar

    % Compute error
    e = d - y; % Shape: scalar

    % Update inverse R matrix
    g_norm_factor = 1 + g' * u;
    Phi_inv = (1 / sqrt(lambda)) * Phi_inv - (g * g') / (sqrt(lambda) * g_norm_factor);

    % Update weights
    w = w + g * e;

    % Update filter state
    filter_state.Phi_inv = Phi_inv;
    filter_state.w = w;
end
