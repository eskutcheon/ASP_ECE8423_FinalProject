% function A = givens_rotation(A, i, M)
%     % Compute and apply a Givens rotation to zero out the (i+1, end) element
%     % Inputs:
%     %   A - Input matrix (augmented)
%     %   i - Row index for the Givens rotation
%     %   M - Number of rows/columns in Phi
%     % Outputs:
%     %   G - Givens rotation matrix (not explicitly used)
%     %   A - Updated matrix after applying the rotation
%     fprintf("A shape: (%d,%d)\n", size(A))
%     % Retrieve the elements to be zeroed and retained
%     a = A(i, i);
%     b = A(i, M+1); % Element in the (i, M+1) position
% 
%     % Compute Givens rotation coefficients
%     r = hypot(a, b); % Compute the hypotenuse directly
%     if r == 0
%         return; % Skip rotation if both a and b are zero
%     end
%     c = a / r; % Cosine component
%     s = -b / r; % Sine component
% 
%     % Apply the rotation to rows i and M+1
%     for col = 1:size(A, 2)
%         temp = c * A(i, col) - s * A(i + 1, col);
%         A(i + 1, col) = s * A(i, col) + c * A(i + 1, col);
%         A(i, col) = temp;
%     end
function [c, s] = givens_rotation(a, b)
    % Compute Givens rotation coefficients to zero out 'b'.
    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            tau = -a / b;
            s = 1 / sqrt(1 + tau^2);
            c = s * tau;
        else
            tau = -b / a;
            c = 1 / sqrt(1 + tau^2);
            s = c * tau;
        end
    end
    % a = A(i, end);   % Element to retain
    % b = A(i+1, end); % Element to zero out
    % 
    % % Compute Givens rotation coefficients
    % r = hypot(a,b);
    % if r == 0
    %     return; % Skip rotation if both a and b are zero
    % end
    % c = a / r;
    % s = -b / r;
    % 
    % % % Create the Givens rotation matrix
    % % G = eye(M + 1);
    % % G([i, i+1], [i, i+1]) = [c, s; -s, c];
    % % fprintf("G shape: (%d,%d)\n", size(G))
    % % % Apply the rotation
    % % A = A * G;
    % % apply rotation to rows i and i+1
    % for col = 1:size(A, 2)
    %     temp = c * A(i, col) - s * A(i+1, col);
    %     A(i+1, col) = s * A(i, col) + c * A(i+1, col);
    %     A(i, col) = temp;
    % end
end
