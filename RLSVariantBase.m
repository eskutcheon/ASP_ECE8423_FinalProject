classdef (Abstract) RLSVariantBase < RLSFilterBase
    properties (Abstract, Access = protected)
        P; % Covariance matrix or its derivative (inverse, square root, etc.)
    end
    
    methods (Abstract)
        % Subclasses must define these
        updateMatrix(obj, input);
        computeGain(obj, input);
    end

    methods
        function adapt(obj, input, desired)
            % General RLS adaptation logic
            gain = obj.computeGain(input);
            aPrioriError = obj.computeAPrioriError(input, desired);
            obj.weights = obj.weights + gain * aPrioriError;
            obj.updateMatrix(input);
        end
    end
end
