% Driver script to compare different adaptive filters

% reproducibility settings
rng default

% parameter initialization
num_samples = 1000;
filter_order = 10;
desired_freq = 0.01;        % frequency of the desired signal
forgetting_factor = 0.8;   % lambda
regularization = 0.01;      % delta
filter_coefs = [1, -0.5];


% synthetic data generation
t = (1:num_samples)';
desired_signal = sin(2*pi*desired_freq*t);      % Clean sinusoidal signal
noise = randn(num_samples, 1);                  % random white Gaussian noise
% Not sure if using MATLAB filter function or a global variable that was
    % never initialized
    % https://www.mathworks.com/help/matlab/ref/filter.html
% essentially x or u - it's the reference noise correlated with the interference but not the desired signal
    % The adaptive filter uses u(n) to estimate e(n) and subtract it from d(n)
reference_noise = filter(filter_coefs, 1, noise);      % Correlated reference noise
% non-stationary reference noise that varies more with time
reference_noise = (1 + 0.5 * sin(2 * pi * 0.01 * t)) .* reference_noise;
noisy_signal = desired_signal + reference_noise;    % Observed noisy signal
%noisy_signal = desired_signal + noise;    % Observed noisy signal

% Compare filters - TODO: remove comments in the actual full test with all
    % filters
%filters = {FilterInvQRDRLS(), LMSFilter()};
%filterNames = {'Inverse QRD-RLS', 'LMS'};

% Initialize RLS filter and noise canceller
rlsFilter = FilterBasicRLS();
%rlsFilter = rlsFilter.initialize(filter_order, forgetting_factor, regularization);
%rlsFilter = FilterKalmanRLS();
%rlsFilter = rlsFilter.initialize(filter_order, forgetting_factor, regularization);


noise_canceller = NoiseCanceller(rlsFilter, filter_order, forgetting_factor, regularization);
% Run the noise canceller
[output, error] = noise_canceller.run(reference_noise, noisy_signal);
%[output, error] = noise_canceller.run(desired_signal, noisy_signal);

% Plot results
figure;
subplot(3, 1, 1);
plot(t, desired_signal, 'b', 'LineWidth', 1.5);
title('Desired Signal');
xlabel('Time (samples)');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(t, noisy_signal, 'r', 'LineWidth', 1.5);
title('Noisy Signal');
xlabel('Time (samples)');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(t, output, 'g', 'LineWidth', 1.5);
title('Filtered Signal');
xlabel('Time (samples)');
ylabel('Amplitude');
%{
subplot(5, 1, 4);
plot(t, reference_noise, 'k', 'LineWidth', 1.5);
title('Noise Reference');
xlabel('Time (samples)');
ylabel('Amplitude');

subplot(5, 1, 5);
plot(t, noise, 'y', 'LineWidth', 1.5);
title('Unfiltered Noise Reference');
xlabel('Time (samples)');
ylabel('Amplitude');
%}
%{
for i = 1:length(filters)
    fprintf('Running experiment with %s filter...\n', filterNames{i});
    canceller = AdaptiveNoiseCanceller(filters{i}, filter_order);
    canceller.runExperiment(noiseSignal, desired_signal);
    canceller.displayMetrics();
end
%}
