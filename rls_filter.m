function [filter_state, y, e] = rls_filter(action, filter_order, N, lambda, delta, varargin)
    % Basic RLS filter implementation
    % FIXME: ensure it's the weighted form - doesn't seem to output anything but the noise at the
    % moment
    persistent w P;
    % NOTE: M = filter length; will use N for the number of samples from here on
    switch action
        case 'initialize'
            % initialize filter weights to 0
            w = zeros(filter_order, 1);
            fprintf("w.shape in rls_filter: (%d,%d)\n", size(w))
            % should be of shape (M, M) according to textbook
            % inverse correlation matrix of the tap-input $u(n)$
                % "may also be viewed as the covariance matrix of RLS estimate
                % $\hat{w}(n)$, normalized with respect to the noise variance $\sigma^2$"
            P = (1 / delta) .* eye(filter_order);
            fprintf("P.shape in rls_filter: (%d,%d)\n", size(P))
            filter_state = struct('w', w, 'P', P, 'lambda', lambda);
        case 'adapt'
            % Adaptation step
            TOL = 1e-3 * delta;
            filter_state = varargin{1};
            x = varargin{2};
            fprintf("x.shape in rls_filter: (%d,%d)\n", size(x))
            d = varargin{3};
            fprintf("d.shape in rls_filter: (%d,%d)\n", size(d))
            M = length(x);
            % gain vector - should be of shape (M, 1) according to textbook
            k = (filter_state.P * x) / (filter_state.lambda + x' * filter_state.P * x + TOL);
            fprintf("k.shape in rls_filter: (%d,%d)\n", size(k))
            % a priori estimation error
            xi = d - filter_state.w' * x;       % a priori estimation error
            % a posteriori weight error (previous error - k*xi)
            %filter_state.eps = filter_state.eps - k * xi';
            fprintf("xi.shape in rls_filter: (%d,%d)\n", size(xi))
            % update weights with a priori error
            filter_state.w = filter_state.w + k * xi';
            fprintf("w.shape in rls_filter: (%d,%d)\n", size(filter_state.w))
            % update inverse covariance matrix
            filter_state.P = (1 / filter_state.lambda) .* (eye(M) - k * x') * filter_state.P;
            fprintf("P.shape in rls_filter: (%d,%d)\n", size(filter_state.P))
            % output signal (filtered signal)
            y = filter_state.w' * x;
            fprintf("y.shape in rls_filter: (%d,%d)\n", size(y))
            % a posteriori estimation error (desired - output)
            e = d - y;
            % seemingly equivalent definition with the conversion factor
            %e = (1 - k' * x) * xi;
            %gamma = e / xi;
            fprintf("e.shape in rls_filter: (%d,%d)\n", size(e))
            filter_state = struct('w', filter_state.w, 'P', filter_state.P, 'lambda', filter_state.lambda);
    end
end