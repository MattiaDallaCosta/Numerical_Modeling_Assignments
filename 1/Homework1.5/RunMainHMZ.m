clear; close all;

%% Path
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%% Data for Test: save the data of the problem.
Data = DataTest('Homework1Point5First');
%Data = DataTest('Homework1Point5Second');

%% Options
Data.visual_graph = 1;      % 1 if you want the plot of u_h and u_ex together
Data.calc_errors = 0;       % 1 if you want to calc many h and compute errors
                            % 0 if you whant to calc for single h and NO errors
                            % graphs
Data.save_sol_images = 1;   % 1 if you whant to save the sol plot

%% number of element definition
nEl = 150;                                       % Single sol calculation
nEl_vec = [1 10 20 30 40 50 60 70 80 80 100 200 1000];   % for many h

%% Main routine
if not (Data.calc_errors)
    [err,sol,fem,D] = MainHMZ(Data,nEl);
end
if (Data.calc_errors)
    
    hVec    = zeros(length(nEl_vec),1);
    eVecL2  = zeros(length(nEl_vec),1);
    eVecH1  = zeros(length(nEl_vec),1);
    
    for i = 1:length(nEl_vec)
    
        nEl = nEl_vec(i);
        fprintf('Running MainHMZ with nEl = %d\n', nEl);
    
        % Run the solver
        [err, sol, fem, D] = MainHMZ(Data, nEl);
    
        % Store errors
        hVec(i)    = fem.h;
        eVecL2(i)  = err.L2;
        eVecH1(i)  = err.H1;
    end
    
    %% Plot Errors
    PlotErr(hVec, eVecL2, eVecH1);
end