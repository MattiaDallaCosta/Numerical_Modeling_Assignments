clear; close all;

%% Path
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing
addpath SemLib


%% Data for Test
%Data = DataTest('Homework2.4');
%Data = DataTest('Homework2.5_1');
Data = DataTest('Homework2.5_2');

% Polynomial degree 
Data.p = 4;

%% Options
Data.surf = 0;          % solution in space-time domain
Data.snapshot = 1;      % movie of the solution 
Data.visual_graph = 0;  % solution at final time T
Data.calc_errors = 0;   % error calculation

% Time steps (bigger for video, to make is shorter)
if (Data.snapshot)
    Data.dt = 0.01;
end

%% No error case
% Input: Data and number of elements
% Output: error, solution, fem and data structure
if not (Data.calc_errors)
    nEl = 10;
    [err,sol,fem,D] = Main(Data,nEl);
end
    
%% Error case
if Data.calc_errors
    nEl = [16 32 64 128];

    hVec   = zeros(1,length(nEl));
    eEtaL2 = zeros(1,length(nEl));
    eEtaH1 = zeros(1,length(nEl));
    eQL2   = zeros(1,length(nEl));
    eQH1   = zeros(1,length(nEl));

    for i = 1:length(nEl)
        [err, sol, fem, D] = Main(Data, nEl(i));
        hVec(i)   = fem.h;
        eEtaL2(i) = err.eta_L2;
        eEtaH1(i) = err.eta_H1;
        eQL2(i)   = err.q_L2;
        eQH1(i)   = err.q_H1;
    end

    figure; clf;

    % L2
    hs = subplot(2,1,1);
    loglog(hVec, hVec.^2, '-+b', 'LineWidth', 2); hold on;
    loglog(hVec, eEtaL2,  '-or', 'LineWidth', 2);
    loglog(hVec, eQL2,    '-sg', 'LineWidth', 2);
    legend(sprintf('h^%i',2), '||\eta-\eta_h||_{L^2}', '||q-q_h||_{L^2}', 'Location', 'best');
    ylabel('L^2-error'); xlabel('h'); hs.FontSize = 12;

    % H1
    hs = subplot(2,1,2);
    loglog(hVec, hVec,    '-+b', 'LineWidth', 2); hold on;
    loglog(hVec, eEtaH1,  '-or', 'LineWidth', 2);
    loglog(hVec, eQH1,    '-sg', 'LineWidth', 2);
    legend(sprintf('h^%i',1), '||\eta-\eta_h||_{H^1}', '||q-q_h||_{H^1}', 'Location', 'best');
    ylabel('H^1-error'); xlabel('h'); hs.FontSize = 12;
end


