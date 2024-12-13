% driver script to compare different adaptive filters based on RLS as part
% of an effort to compare against a final implementation of QRD-RLS

function filter_state = initialize_filter_state(filter_func, M, N, lambda, delta, x, d)
    switch func2str(filter_func)
        case 'rls_filter'
            % RLS Filter state
            filter_state = struct( ...
                'w', zeros(M, 1), ...               % filter weights
                'P', (1 / delta) .* eye(M), ...     % M x M inverse covariance matrix
                'lambda', lambda ...                % forgetting factor
            );

        case 'kalman_rls_filter'
            x0_white = x(M:-1:1) - mean(x(M:-1:1));
            % Kalman Filter state
            filter_state = struct( ...
                'x_hat', x(M:-1:1), ...       % initial state
                'K', mean(x0_white * x0_white') .* eye(M), ...   % initial error covariance matrix
                'Q1', var(x) * eye(M), ...             % process noise covariance
                'Q2', var(d) * eye(N), ...             % measurement noise covariance
                'F', (1 - delta) .* eye(M), ...          % state transition matrix with weight decay
                'C', ceil(sprand(N, M, 0.5)) ...   % measurement equation matrix
            ) % 'C', repmat(x, 1, M) ...

        otherwise
            error('Unknown filter function');
    end
end


function [output_signal, error_signal] = run_noise_canceller(filter_func, M, lambda, delta, ref_noise, noisy_signal)
    % initialize the output array and error array
    N = length(noisy_signal);
    output_signal = zeros(N, 1);
    error_signal = zeros(N, 1);
    fprintf("output and error shapes in noise canceller: (%d,%d)\n", size(output_signal))
    % initialize filter state variables as a struct
    %filter_state = filter_func('initialize', M, N, lambda, delta, ref_noise(M:-1:1));
    %filter_state = filter_func('initialize', M, lambda, delta, ref_noise);
    filter_state = initialize_filter_state(filter_func, M, N, lambda, delta, ref_noise, noisy_signal);
    % adapt filter by iterating over samples
    for iter = 1:10
        switch func2str(filter_func)
            case 'rls_filter'
                for n = M:N
                    %input = ref_noise(n:-1:(n - M + 1));
                    input_signal = ref_noise(n:-1:n-M+1);
                    fprintf("input shape in noise canceller: (%d,%d)\n", size(input_signal))
                    desired = noisy_signal(n);
                    fprintf("desired signal shape in noise canceller: (%d,%d)\n", size(desired))
                    [filter_state, output_signal(n), error_signal(n)] = filter_func(filter_state, input_signal, desired);
                    % TODO: after implenting new logger, call it here
                end
            case 'kalman_rls_filter'
                [filter_state, output_signal, error_signal] = filter_func(filter_state, ref_noise, noisy_signal);
            otherwise
                error('Unknown filter function');
        end
    end
end


function plot_results(t, desired_signal, noisy_signal, output_signal, filter_order, filter_name)
    % Plot results
    figure;
    subplot(3, 1, 1);
    t1 = t(filter_order:end);
    % plot of underlying target signal buried by noise
    plot(t1, desired_signal(filter_order:end), 'b', 'LineWidth', 1.5);
    title(['Desired Signal - ', filter_name]);
    xlabel('Time (samples)');
    ylabel('Amplitude');
    % plot of signal corrupted by noise
    subplot(3, 1, 2);
    plot(t1, noisy_signal(filter_order:end), 'r', 'LineWidth', 1.5);
    title('Noisy Signal');
    xlabel('Time (samples)');
    ylabel('Amplitude');
    % filtered output signal plot
    subplot(3, 1, 3);
    plot(t1, output_signal(filter_order:end), 'g', 'LineWidth', 1.5);
    title('Filtered Signal');
    xlabel('Time (samples)');
    ylabel('Amplitude');
end



% reproducibility settings
rng default

% parameter initialization
num_samples = 1000;
% TODO: figure out if it would be accurate to call this the number of tap weights like the textbook
filter_order = 10;
desired_freq = 0.01;        % frequency of the desired signal
forgetting_factor = 0.99;    % lambda
regularization = 1e-3;      % delta
filter_coefs = [1, -0.5];
%filter_coefs = 2*rand(filter_order, 1) - 1 % M-length vector in [-1,1]

% synthetic data generation
t = (1:num_samples)';
desired_signal = sin(2*pi*desired_freq .* t);    % Clean sinusoidal signal
noise = randn(num_samples, 1);      % random white Gaussian noise
% essentially x or u - it's the reference noise correlated with the interference but not the desired signal
    % The adaptive filter uses u(n) to estimate e(n) and subtract it from d(n)
    % https://www.mathworks.com/help/matlab/ref/filter.html
ref_noise = filter(filter_coefs, 1, noise);      % Correlated reference noise
% non-stationary reference noise that varies more with time
ref_noise = (1 + 0.5 .* sin(2 * pi * desired_freq .* t)) .* ref_noise;
noisy_signal = desired_signal + ref_noise;    % Observed noisy signal
%noisy_signal = desired_signal + noise;    % Observed noisy signal
% Run experiments with different filters
%filters = {@rls_filter, @kalman_rls_filter};
%filter_names = {"Basic RLS", "Kalman RLS"};
filters = {@kalman_rls_filter};
filter_names = {"Kalman RLS"};

for i = 1:length(filters)
    fprintf("Running %s...\n", filter_names{i});
    [output, err] = run_noise_canceller(filters{i}, filter_order, forgetting_factor, regularization, ref_noise, noisy_signal);
    output_signal = output;
    plot_results(t, desired_signal, ...
                    noisy_signal, ...
                    output_signal, filter_order, filter_names{i});
    %return;
end
