function K = lookupHalfSaturationConstants(obj, components)
    arguments
        obj MPC.ecological.Reaction
        components MPC.ecological.Component
    end
    names = [components.Name];

    numReactions = length(obj);
    numComponents = length(names);
    
    K = zeros(numComponents, numReactions);
    for i = 1:length(obj)
        K(:, i) = lookup(obj(i).HalfSaturationConstants, names, FallbackValue=nan);
    end
end