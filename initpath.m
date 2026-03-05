function initpath()
    if isfolder("./src")
        addpath(genpath("./src"))
    end
end