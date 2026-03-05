classdef Model
    properties
        Components (:, 1)
        
        BiofilmPorosity
        OsmosisRate
        DetachmentFunction

        Reactions = MPC.ecological.Reaction.empty;
        CohesionSubModel = MPC.cohesion.CahnHilliardModel()
    end
    
    properties (Dependent, Hidden)
        Particles
        Liquids
        
        StoichiometricCoefficients
        HalfSaturationConstants
        Order

        StoichiometricMatrixParticles
        StoichiometricMatrixLiquids
    end

    methods
        function obj = Model(input)
            arguments
                input.Components = MPC.ecological.Component.empty;
                input.Kappa double = [];
                input.Zeta0 double = [];
                input.Zeta1 double = [];
                input.MobilityFunction function_handle = @(u) u .* (1 - u);
                input.PotentialGradient function_handle = @(u) 0.25*((u.^2).*((1 - u).^2));

                input.DetachmentFunction = @(v, qnom) sqrt(abs(v)/qnom);
                input.OsmosisRate = 1e-5;
                input.BiofilmPorosity = .99;

                input.Reactions = MPC.ecological.Reaction.empty;               
                input.Preset = string.empty;
            end
            switch input.Preset
                case string.empty
                    propNames = string(properties(obj));
                    for field = propNames(~contains(propNames, "SubModel")).'
                        obj.(field) = input.(field);
                    end
                    for superField = propNames(contains(propNames, "SubModel")).'
                        for field = string(fieldnames(obj.(superField)).')
                            inputValue = input.(field);
                            if ~isempty(inputValue)
                                obj.(superField).(field) = inputValue;
                            end
                        end
                    end
                case {"Lund", "Rosenqvist"}
                    obj = MPC.presets.modelLund();
                otherwise
                    error("Unknown model preset.");
            end
        end
        
        function p = get.Particles(obj)
            p = obj.Components(arrayfun(@(C)isa(C, 'MPC.ecological.Particle'), obj.Components));
        end

        function l = get.Liquids(obj)
            l = obj.Components(arrayfun(@(C)isa(C, 'MPC.ecological.Liquid'), obj.Components));
        end

        function sigma = get.StoichiometricCoefficients(obj)
            sigma = obj.Reactions.lookupStoichiometricCoefficients(obj.Components);
        end
        
        function sigma = get.StoichiometricMatrixParticles(obj)
            sigma = obj.Reactions.lookupStoichiometricCoefficients(obj.Particles);
        end

        function sigma = get.StoichiometricMatrixLiquids(obj)
            sigma = obj.Reactions.lookupStoichiometricCoefficients(obj.Liquids);
        end

        function K = get.HalfSaturationConstants(obj)
            K = obj.Reactions.lookupHalfSaturationConstants(obj.Components);
        end

        function p = get.Order(obj)
            p = obj.Reactions.lookupOrder(obj.Components);
        end
        
        mu = computeReactionRates(obj, temperature);
    end
end