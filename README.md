# AMMONIA
multiphAse continuuM Model of slOw saNd fIltrAtion

## Example

```matlab
clc; clear; restoredefaultpath;
initpath;
import MPC.*;

filter = addGridPoints(SandFilter(), 200);
model = load("./data/LundMPCModel.mat").model;
inflow = [1e-4, 2e-4, 0, 0, 1e-5, 1e-6, 1e-9, 0, 1e-5];

filterState = State(filter, model);
results = filterState.simulate(InflowConcentrations=inflow, SimulationTime=10.0);
```
