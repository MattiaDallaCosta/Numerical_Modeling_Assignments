clear; close all;

%% Path
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%% Data for Test: save the data of the problem.
Data = DataTest('Homework1Point2');

%% Options
Data.visual_graph = 0; % 1 if you want the plot of u_h and u_ex together
Data.calc_errors = 1;  % 1 if you want to calc many for h and compute errors
                       % 0 if you whant to calc for single h and NO errors
                       % graphs
Data.save_sol_images = 1;  % 1 if you whant to save the sol plot

%% number of element definition
nEl = 10;                                       % Single sol calculation

nEl_min   = 1;         % coarsest mesh
nEl_max   = 1000;     % finest mesh
n_levels  = 200;      % number of elements

% Logarithmically spaced nEL:
nEl_vec = unique( round( logspace(log10(nEl_min), log10(nEl_max), n_levels) ) );

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