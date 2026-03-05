classdef Particle < MPC.ecological.Component & matlab.mixin.CustomDisplay
    properties
        AttachmentRate double = 0.0 % [1/day]
        Attenuation double = 0.0 % [m^2/kg]
    end
    methods
        function obj = Particle(input)
            arguments
                input.Name string = string.empty
                input.Density double = []
                input.Dispersivity double = []
                input.AttachmentRate double = []
                input.TransportRate double = []
                input.Attenuation double = []
            end
            extra_fields = ["AttachmentRate", "Attenuation"];
            input_args = namedargs2cell(rmfield(input, extra_fields));
            obj = obj@MPC.ecological.Component(input_args{:});
            for field = extra_fields;
                input_value = input.(field);
                if ~isempty(input_value)
                    obj.(field) = input_value;
                end
            end
        end
    end

    methods (Access = protected)
        function displayNonScalarObject(objArray)
            dimStr = matlab.mixin.CustomDisplay.convertDimensionsToString(objArray);
            cName = matlab.mixin.CustomDisplay.getClassNameForHeader(objArray);
            MPC.ecological.Component.displayHomogeneousNonScalarObject(objArray, dimStr, cName);
        end
    end
end