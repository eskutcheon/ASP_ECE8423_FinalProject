function w = back_substitution(Phi, p)
    % Solve the triangular system Phi' * w = p
    % Inputs:
    %   Phi - Upper triangular Cholesky factor (M x M)
    %   p   - Projection vector (M x 1)
    % Output:
    %   w   - Weight vector (M x 1)
    M = length(p);
    w = zeros(M, 1);
    % Solve using back-substitution
    for i = M:-1:1
        w(i) = (p(i) - Phi(i, i+1:end) * w(i+1:end)) / Phi(i, i);
    end
end