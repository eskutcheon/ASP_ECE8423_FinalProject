function [filter_state, y, e] = rls_filter(filter_state, x, d)
    % Basic RLS filter implementation
    % NOTE: M = filter length; will use N for the number of samples from here on
    TOL = 1e-8;
    %%%fprintf("x.shape in rls_filter: (%d,%d)\n", size(x))
    %%%fprintf("d.shape in rls_filter: (%d,%d)\n", size(d))
    M = length(x);

    % output signal (filtered signal)
    y = filter_state.w' * x;
    %%%fprintf("y.shape in rls_filter: (%d,%d)\n", size(y))
    
    % a priori estimation error
    xi = d - y;       % a priori estimation error
    %%%fprintf("xi.shape in rls_filter: (%d,%d)\n", size(xi))

    % gain vector - should be of shape (M, 1) according to textbook
    k = (filter_state.P * x) / (filter_state.lambda + x' * filter_state.P * x + TOL);
    %%%fprintf("k.shape in rls_filter: (%d,%d)\n", size(k))

    % update weights with a priori error
    filter_state.w = filter_state.w + k * xi';
    %%%fprintf("w.shape in rls_filter: (%d,%d)\n", size(filter_state.w))
    % update inverse covariance matrix
    filter_state.P = (1 / filter_state.lambda) .* (eye(M) - k * x') * filter_state.P;
    %%%fprintf("P.shape in rls_filter: (%d,%d)\n", size(filter_state.P))

    % seemingly equivalent definition with the conversion factor
    %e = (1 - k' * x) * xi;
    %gamma = e / xi;
    % a posteriori weight error (previous error - k*xi)
    y = filter_state.w' * x;
    e = d - y;
    % NOTE: still not sure if this is necessary, since the actual
    % parameters that need to be updated are done through the struct anyway
    %filter_state = struct('w', filter_state.w, 'P', filter_state.P, 'lambda', filter_state.lambda);
end