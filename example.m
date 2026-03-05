clc; clear; restoredefaultpath; 
initpath; 
import MPC.*;

filter = addGridPoints(SandFilter(), 200);
model = load("./LundMPCModel.mat").model;
state = State(filter, model);

tic
results = simulate(state, ...
    InflowConcentrations=[1e-3, 1e-3, 0, 1e-5, 1e-2, 1e-2, 1e-5, 0, 1e-4] / 10, ...
    SimulationTime=1.0, TimeStep=1e-6);
toc