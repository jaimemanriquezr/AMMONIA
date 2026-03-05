function mu = computeReactionRates(obj, temperature)
    arguments
        obj MPC.Model
        temperature double
    end
    mu = arrayfun(@(r) r.computeRate(temperature), obj.Reactions);
end