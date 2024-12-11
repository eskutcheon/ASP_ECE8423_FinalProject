classdef NoiseCanceller
    properties
        filter;         % Instance of some adaptive filter class
        logger;         % MetricsLogger instance
        filter_order;    % Number of filter tap-weights
    end
    
    methods
        function obj = NoiseCanceller(filter, filter_order, varargin)
            % Constructor - sets the adaptive filter and filter order
            obj.filter = filter;
            obj.filter_order = filter_order;
            % Initialize filter with variable arguments
            obj.filter = obj.filter.initialize(filter_order, varargin{:});
            % obj.logger = MetricsLogger();
        end
        
        
        function [output_signal, filter_error] = run(obj, noise_reference, noisy_signal)
            % Initialize the filter, output array, and error array
            M = obj.filter_order;
            % TODO: add safeguards on the input arguments
            num_samples = length(noisy_signal);
            output_signal = zeros(num_samples, 1);
            filter_error = zeros(num_samples, 1);
            % Adapt filter over the input signal by iterating over samples
            for i = M:num_samples
                input = noise_reference(i:-1:(i - M + 1));
                %input = noisy_signal(i:-1:(i - M + 1));
                desired = noisy_signal(i);
                [obj.filter, output_signal(i), filter_error(i)] = obj.filter.adapt(input, desired);
                %output_signal(i) = output_signal(i) - filter_error(i);
                %fprintf("[output_signal(i), filter_error(i)] = (%f, %f)\n", [output_signal(i), filter_error(i)])
            end
        end
        % TODO: implement new logger as a struct of cells later
    end
end