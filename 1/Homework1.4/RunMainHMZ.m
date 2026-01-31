clear; close all;

% Path
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


% Data for Test: save the data of the problem.
Data = DataTest('Homework1Point2');

% Options
Data.visual_graph = 0; % 1 if you want the plot of u_h and u_ex together
Data.calc_errors = 1;  % 1 if you want to calc many for h and compute errors
                       % 0 if you whant to calc for single h and NO errors
                       % graphs
Data.save_sol_images = 1;  % 1 if you whant to save the sol plot
Data.plot_matrixes = 1;

% Number of modes to consider
nmodes = 5;

% number of elements on which calculate the different dispersion errors
%nEl_min   = 1;        % coarsest mesh
%nEl_max   = 1000;     % finest mesh
%n_levels  = 100;      % number of elements
%nEl_vec = unique( round( logspace(log10(nEl_min), log10(nEl_max), n_levels) ) );
nEl_vec = [5 7 10 25 50 100 200 400 800];

% Compute analytical k for each one of the modes
n = (1:nmodes);
k_exact = zeros(nmodes,1);          % Store analytical wavenumbers
k_exact(:,1) = (2*n - 1)*pi/(2);    % Compute analytical wavenumbers

% prepare arrays to store dispertion ratios and relative h/lambda
dispersion_ratio = zeros(nmodes*length(nEl_vec),1);
h_over_lambda  = zeros(nmodes*length(nEl_vec),1);

for i = 1:length(nEl_vec)

    nEl = nEl_vec(i);
    fprintf('Running MainHMZ with nEl = %d\n', nEl);

    % Run the solver: it will compute the values k/k_h and h/lambda for all
    % the nmodes selected. The nmodes values caluclated are then saved in
    % the two arrays previusly created
    [dispersion_ratio(nmodes*(i-1)+1:nmodes*(i)), h_over_lambda(nmodes*(i-1)+1:nmodes*(i))] = MainHMZ(Data, nEl, k_exact, nmodes);
end
    
%% Plot Errors
PlotDisp(dispersion_ratio, h_over_lambda, nmodes, Data);