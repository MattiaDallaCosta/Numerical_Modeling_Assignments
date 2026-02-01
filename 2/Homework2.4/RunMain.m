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
Data = DataTest('Homework2.4');

%% Options
Data.surf = 0;          % solution in space-time domain
Data.snapshot = 0;      % movie of the solution 
Data.visual_graph = 1;  % solution at final time T
Data.calc_errors = 0;   % error calculation

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
    eVecL2 = zeros(1,length(nEl));
    eVecH1 = zeros(1,length(nEl));

    for i = 1:length(nEL)
        [err,sol,fem,D] = Main(Data,nEl(i));
        hVec(1,i)   = fem.h;    % Mesh size
        eVecL2(1,i) = err.L2;   % L2 norm error
        eVecH1(1,i) = err1.H1;  % H1 norm error
    end
    
    %Plot
    
    hs = subplot(2,1,1);
    loglog(hVec,hVec.^2,'-+b','Linewidth',2); hold on;
    loglog(hVec,eVecL2,'-or','Linewidth',2);
    legend(sprintf('h^%i',2),'||u-u_h||_{L^2}');
    ylabel('L^2-error');
    xlabel('h');
    hs.FontSize = 12;
    
    hs = subplot(2,1,2);
    loglog(hVec,hVec,'-+b','Linewidth',2); hold on;
    loglog(hVec,eVecH1,'-or','Linewidth',2);
    legend(sprintf('h^%i',1),'||u-u_h||_{H^1}');
    ylabel('H^1-error')
    xlabel('h');
    hs.FontSize = 12;
end

