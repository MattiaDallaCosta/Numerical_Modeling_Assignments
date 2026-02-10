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
%Data = DataTest('Homework2.5_1');
%Data = DataTest('Homework2.5_2');

% Polynomial degree 
Data.p = 2;

%% Options
Data.surf = 0;            % solution in space-time domain (3d)
Data.snapshot = 0;        % movie of the solution
Data.visual_graph = 0;    % solution at final time T
Data.calc_errors = 1;     % error calculation
Data.save_sol_images = 1; % Save image

% Skip some steps in the video to make it go faster
if (Data.snapshot)
    Data.frameSkip = 100;
end

% The following is both the time at which shows the sol and computes the
% error
Data.PicTime = 0.5;

%% No error case
% Input: Data and number of elements
% Output: error, solution, fem and data structure
if not (Data.calc_errors)
    nEl = 10;
    [err,sol,fem,D] = Main(Data,nEl);
end
    
%% Error case
if Data.calc_errors
    %nEl = [5 6 7 8 9 10 12 15 20 22 25 50 100 200];
    nEl = [5 6 7 8 9 10 12 15 17 20 25 50 100 200];
    %nEl = [5 10 15 50];

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
    
    s = diff(log(eEtaL2)) ./ diff(log(hVec));    % with h ascending
    pSpace = Data.p + 1;

    % Use slopes only away from extremes (optional)
    I = 2:length(s)-1;                           % avoid endpoints
    [~, j] = min(abs(s(I) - pSpace));
    kAnchor = I(j) + 1;                          % slope s(k-1) relates to point k
    
    C_L2  = eEtaL2(kAnchor) / (hVec(kAnchor)^pSpace);
    refL2 = C_L2 * hVec.^pSpace;

    % Build reference lines scaled to match last error point
    %C_L2 = eEtaL2(end) / (hVec(end)^(Data.p+1));
    %refL2 = C_L2 * hVec.^(Data.p+1);
    
    % % Compute rate
    % %---------------------- L2 error ----------------------%
    % loglog(hVec, refL2, '-+b', 'LineWidth', 2); hold on;
    % loglog(hVec, eVecL2, '-or', 'LineWidth', 2);
    % pfitL2 = polyfit(log(hVec), log(eVecL2), 1);
    % rateL2 = pfitL2(1);

    fs         = 35;
    lw         = 1.9;
    fnt        = 'Times New Roman';
    interp_tex = 'latex';
    
    clr = [ ...
        0.00  0.45  0.74;   %  1  blue
        0.85  0.33  0.10;   %  2  orange
        0.47  0.67  0.19;   %  3  green
        0.49  0.18  0.56;   %  4  purple
        0.93  0.69  0.13;   %  5  yellow
        0.30  0.75  0.93;   %  6  cyan
        0.64  0.08  0.18;   %  7  dark red
        0.25  0.25  0.25;   %  8  dark gray
        0.00  0.60  0.50;   %  9  teal
        0.75  0.00  0.75;   % 10  magenta
        0.60  0.60  0.60];  % 11  light gray


    fig3 = figure('Units','centimeters','Position',[2 2 30 16], ...
              'Color','w','PaperPositionMode','auto');
    hold on; box on; grid on;
    set(gca, 'GridLineStyle','-', 'GridAlpha',0.15);
    
    loglog(hVec, refL2,   '-+', 'LineWidth',lw, 'MarkerSize',8, 'Color',clr(1,:));
    loglog(hVec, eEtaL2,  '-o', 'LineWidth',lw, 'MarkerSize',6, 'Color',clr(2,:));
    set(gca,'XScale','log')
    set(gca,'YScale','log')

    legend({sprintf('$C h^{%d}$', Data.p+1), ...
            '$\|\eta-\eta_h\|_{L^2}$'}, ...
           'Interpreter',interp_tex, ...
           'Location','best', ...
           'FontName',fnt, ...
           'FontSize',fs-8, ...
           'Box','on');
    
    xlabel('$h$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    ylabel('$L^2$-error', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    title(sprintf('$L^2$ errors at time $%.1f$\\,s with $p = %d$ and $\\Delta t = %.2f$\\,ms', ...
          Data.PicTime, Data.p, Data.dt*1000), ...
          'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    
    set(gca, 'FontName',fnt, 'FontSize',fs-4, 'TickLabelInterpreter','latex');
    
    if Data.save_sol_images
        filename = sprintf('ErrorPlot_q_%s_p_%g_dt_%g_t_%.1f.pdf', ...
                           Data.name, Data.p, Data.dt, Data.PicTime);
        exportgraphics(fig3, fullfile('Plots', filename), 'ContentType','vector');
    end
    % figure; clf;
    % 
    % % L2
    % % hs = subplot(2,1,1);
    % loglog(hVec, refL2, '-+b', 'LineWidth', 1.8); hold on;
    % loglog(hVec, eEtaL2,  '-or', 'LineWidth', 1.8);
    % %loglog(hVec, eQL2,    '-sg', 'LineWidth', 1.8);
    % %legend(sprintf('h^%i',2), '||\eta-\eta_h||_{L^2}', '||q-q_h||_{L^2}', 'Location', 'best', FontSize=16);
    % legend(sprintf('Ch^%i',Data.p+1), '||\eta-\eta_h||_{L^2}', 'Location', 'best', FontSize=16);
    % ylabel('L^2-error', FontSize=17); xlabel('h', FontSize=17);
    % title(sprintf('L2 errors at time %.1f s with p = %.0f and dt = %.2f ms', Data.PicTime, Data.p, Data.dt*1000), 'FontSize', 16);
    % 
    % % H1
    % % hs = subplot(2,1,2);
    % % loglog(hVec, hVec,    '-+b', 'LineWidth', 1.8); hold on;
    % % loglog(hVec, eEtaH1,  '-or', 'LineWidth', 1.8);
    % % loglog(hVec, eQH1,    '-sg', 'LineWidth', 1.8);
    % % legend(sprintf('h^%i',1), '||\eta-\eta_h||_{H^1}', '||q-q_h||_{H^1}', 'Location', 'best', FontSize=16);
    % % ylabel('H^1-error', FontSize=17); xlabel('h', FontSize=17);
    % % title(sprintf('H^1 errors at time %.1f s with p = %.2f and h = %.2f', Data.PicTime, Data.p, 1/fem.ne), 'FontSize', 17);
    % 
    % 
    % if (Data.save_sol_images)
    %     filename = sprintf('ErrorPlot_q_%s_p_%g_dt_%g_t_%.1f.png', Data.name, Data.p, Data.dt, Data.PicTime);
    %     exportgraphics(gcf, fullfile('Plots', filename), 'Resolution', 300);
    % end
end


