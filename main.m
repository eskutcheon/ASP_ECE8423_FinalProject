% driver script to compare different adaptive filters based on RLS as part
% of an effort to compare against a final implementation of QRD-RLS

function [output_signal, error_signal] = run_noise_canceller(filter_func, M, lambda, delta, ref_noise, noisy_signal)
    % initialize the output array and error array
    N = length(noisy_signal);
    output_signal = zeros(N, 1);
    error_signal = zeros(N, 1);
    fprintf("output and error shapes in noise canceller: (%d,%d)\n", size(output_signal))
    % initialize filter state variables as a struct
    filter_state = filter_func('initialize', M, N, lambda, delta, ref_noise(M:-1:1));
    %filter_state = filter_func('initialize', M, lambda, delta, ref_noise);
    % adapt filter by iterating over samples
    for i = M:N
        input = ref_noise(i:-1:(i - M + 1));
        %input = ref_noise;
        fprintf("input shape in noise canceller: (%d,%d)\n", size(input))
        desired = noisy_signal(i);
        fprintf("desired signal shape in noise canceller: (%d,%d)\n", size(desired))
        [filter_state, output_signal(i), error_signal(i)] = filter_func('adapt', M, N, lambda, delta, filter_state, input, desired);
        % TODO: after implenting new logger, call it here
    end
end


function plot_results(t, desired_signal, noisy_signal, output, filter_name)
    % Plot results
    figure;
    subplot(3, 1, 1);
    % plot of underlying target signal buried by noise
    plot(t, desired_signal, 'b', 'LineWidth', 1.5);
    title(['Desired Signal - ', filter_name]);
    xlabel('Time (samples)');
    ylabel('Amplitude');
    % plot of signal corrupted by noise
    subplot(3, 1, 2);
    plot(t, noisy_signal, 'r', 'LineWidth', 1.5);
    title('Noisy Signal');
    xlabel('Time (samples)');
    ylabel('Amplitude');
    % filtered output signal plot
    subplot(3, 1, 3);
    plot(t, output, 'g', 'LineWidth', 1.5);
    title('Filtered Signal');
    xlabel('Time (samples)');
    ylabel('Amplitude');
end



% reproducibility settings
rng default

% parameter initialization
num_samples = 1000;
% TODO: figure out if it would be accurate to call this the number of tap weights like the textbook
filter_order = 20;
desired_freq = 0.01;        % frequency of the desired signal
forgetting_factor = 0.95;    % lambda
regularization = 1e-3;      % delta
filter_coefs = [1, -0.5];
%filter_coefs = 2*rand(filter_order, 1) - 1 % M-length vector in [-1,1]

% synthetic data generation
t = (1:num_samples)';
desired_signal = sin(2*pi*desired_freq .* t);    % Clean sinusoidal signal
noise = randn(num_samples, 1);                % random white Gaussian noise
% essentially x or u - it's the reference noise correlated with the interference but not the desired signal
    % The adaptive filter uses u(n) to estimate e(n) and subtract it from d(n)
    % https://www.mathworks.com/help/matlab/ref/filter.html
ref_noise = filter(filter_coefs, 1, noise);      % Correlated reference noise
% non-stationary reference noise that varies more with time
ref_noise = (1 + 0.5 .* sin(2 * pi * desired_freq .* t)) .* ref_noise;
noisy_signal = desired_signal + ref_noise;    % Observed noisy signal
%noisy_signal = desired_signal + noise;    % Observed noisy signal
% Run experiments with different filters
filters = {@rls_filter, @kalman_rls_filter};
filter_names = {"Basic RLS", "Kalman RLS"};

for i = 1:length(filters)
    fprintf("Running %s...\n", filter_names{i});
    [output, error] = run_noise_canceller(filters{i}, filter_order, forgetting_factor, regularization, ref_noise, noisy_signal);
    plot_results(t, desired_signal, noisy_signal, output, filter_names{i});
end
