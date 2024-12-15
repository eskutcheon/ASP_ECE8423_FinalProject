function [filter_state, y, e] = qrd_rls_filter(filter_state, u, d)
    % Base QRD-RLS Algorithm
    %fprintf("shape of A: (%d)\n", length(A))
    M = size(u, 1);
    Phi = filter_state.Phi; % upper triangular Cholesky factor
    w = filter_state.w; % filter weights
    p = filter_state.p;
    TOL = 1e-8;
    lambda_sqrt = sqrt(filter_state.lambda);
    %fprintf("Phi (initial) shape: (%d,%d)\n", size(Phi))
    %fprintf("u shape: (%d,%d)\n", size(u))
    % rotate with Givens rotations to get updates for Phi, p, and w
    prearray = [
        lambda_sqrt * Phi, u;
        lambda_sqrt * p', d;
        zeros(1,M), 1
    ];
    %fprintf("prearray shape: (%d,%d)\n", size(prearray))
    %fprintf("shape of R: %d\n", length(R))
    % Update R using Givens rotations to incorporate input_signal
    % for i = 1:M
    %     % Generate Givens rotation matrix to zero out the lower elements
    %     prearray = givens_rotation(prearray, i, M);
    % end
    % Apply Givens rotations to zero out lower triangular elements
    for i = 1:M
        % Compute Givens rotation coefficients to zero out prearray(i+1:end, i)
        %[c, s] = givens_rotation(prearray(i, i), prearray(i+1, i));
        [c, s] = givens_rotation(prearray(i, M+1), prearray(M+1, M+1));
        % Apply Givens rotation to the relevant rows
        G = [c, s; -s, c]; % 2x2 Givens rotation matrix
        prearray([i, M+1], :) = G * prearray([i, M+1], :); % Rotate rows i and i+1
    end
    disp(prearray)
    if sum(isnan(prearray)) > 0
        fprintf("nans detected in prearray: %d\n", sum(isnan(prearray)))
        error("NaNs detected! Terminating...")
    end
    %disp(prearray)
    % Extract the updated Phi, p, and w from the rotated array
    filter_state.Phi = prearray(1:M, 1:M);
    filter_state.p = prearray(M+1, 1:M)';
    gamma_sqrt = prearray(M+1, M+1); % Bottom-right element
    %filter_state.w = inv(filter_state.Phi + TOL * eye(M))' * filter_state.p;
    filter_state.w = back_substitution(filter_state.Phi, filter_state.p);
    %disp(filter_state.w)
    % %fprintf("shape of w: (%d,%d)\n", size(w))
    if sum(isnan(filter_state.w))
        fprintf("nans detected in w: %d\n", sum(isnan(filter_state.w)))
        error("NaNs detected! Terminating...")
    end
    % Filter output and error
    y = filter_state.w' * u; % Filtered output
    e = d - y; % Error signal
end
