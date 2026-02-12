% =========================================================================
% FINITE VOLUME - 1D SHALLOW WATER EQUATIONS
% FIRST ORDER GODUNOV + SECOND ORDER (Lax-Wendroff + MC limiter)
% =========================================================================
%
% We solve the 1D SWE system:
%
%   U_t + F(U)_x = 0
% with
%   - U = [h; q]   with q = h*u
%   - F(U) = [ q ;  q^2/h + (1/2)*g*h^2 ]

clearvars; clc; close all;

%% Outputs
OverlayPlots = 0; % Plot overlay
HeatmapPlots = 0;
ConvergencePlots = 0;
VideoSol = 1;

%% GLOBAL PARAMETERS
g   = 9.81;          % Gravitational acceleration [m/s^2]
L   = 1;             % Domain length [m]
N   = 100;           % Number of cells
%CFL = 0.45;          % CFL number (< 0.5 for stability)
CFL = 0.2;
T   = 1;          % Final time [s]

% Spatial discretization: 
Delta_x = L / N;                                % dx = L/N
x_edges = linspace(0, L, N+1);                  % edges of each cell
xc = (x_edges(1:end-1) + x_edges(2:end)) / 2;   % centers of each cell

% Formatting constants
fs   = 22;           % Font size
lw   = 1.9;          % Line width
fnt  = 'Times New Roman';
interp_tex = 'latex';

% Definition of I.C. U(x,0). q(x,0) = 0 is zero in all cases
cases = struct();
cases(1).tag = 'a'; cases(1).name = 'Dam Break';
cases(1).hL = 1.0;  cases(1).hR = 0.1;

cases(2).tag = 'b'; cases(2).name = 'Reverse Dam Break';
cases(2).hL = 0.1;  cases(2).hR = 1.0;

cases(3).tag = 'c'; cases(3).name = 'Hydraulic Jump';
cases(3).hL = 1.0;  cases(3).hR = 0.5;

% Snapshot times for time-evolution plots
snap_times = [0.00, 0.05, 0.10, 0.20, 0.50, 0.60, 1];

% IC ghost color
ic_color = [0.65 0.65 0.65];

% Create output folder
if ~exist('Plots','dir'), mkdir('Plots'); end

%% RUN BOTH SCHEMES FOR ALL CASES
res_1st = cell(1,3);
res_2nd = cell(1,3);

for tc = 1:3
    fprintf('=== Case (%s): %s ===\n', cases(tc).tag, cases(tc).name);

    h0 = zeros(1,N);
    h0(xc < 0.5) = cases(tc).hL;
    h0(xc >= 0.5) = cases(tc).hR;
    q0 = zeros(1,N);

    % First order Godunov
    fprintf('  1st order ... ');
    [hs1, qs1, ts1, ha1, qa1, ta1] = Godunov1st(h0, q0, xc, Delta_x, g, CFL, T, snap_times);
    res_1st{tc} = struct('hs',hs1,'qs',qs1,'ts',ts1,'ha',ha1,'qa',qa1,'ta',ta1,'h0',h0,'q0',q0);
    fprintf('done (%d steps)\n', length(ta1));

    % Second order Godunov (Lax-Wendroff)
    fprintf('  2nd order (LW+MC) ... ');
    [hs2, qs2, ts2, ha2, qa2, ta2] = Godunov2nd(h0, q0, xc, Delta_x, g, CFL, T, snap_times);
    res_2nd{tc} = struct('hs',hs2,'qs',qs2,'ts',ts2,'ha',ha2,'qa',qa2,'ta',ta2,'h0',h0,'q0',q0);
    fprintf('done (%d steps)\n', length(ta2));
end

%% OVERLAY PLOTS
% For each case: 1 figure, 2 subplots (h top, q bottom).
% All three time instants overlaid; each instant has its own colour,
% 1st order = solid, 2nd order = dashed + sparse marker.
t_overlay = [0.10, 0.20, 1.00];

% --- Per-instant styles (colour / marker) ---
clr_overlay = [ ...
    0.00  0.45  0.74;    % blue   – t1
    0.85  0.33  0.10;    % orange – t2
    0.47  0.67  0.19];   % green  – t3
mrk_overlay = {'o', 's', 'd'};          % marker per instant
mrk_every   = max(1, round(N / 20));    % plot a marker every ~20 pts

if OverlayPlots
    for tc = 1:3
        fig = figure('Units','centimeters','Position',[2 2 30 22], ...
                     'Color','w','PaperPositionMode','auto');

        h0_ic = res_1st{tc}.h0;
        q0_ic = res_1st{tc}.q0;
        dt_avg = mean(diff(res_1st{tc}.ta));

        % ---- subplot 1: h(x,t) ----
        subplot(2,1,1); hold on;
        plot(xc, h0_ic, '-', 'LineWidth', 1.2, 'Color', ic_color);   % IC ghost
        leg_h = {'IC ($t=0$)'};

        for jj = 1:length(t_overlay)
            t_now = t_overlay(jj);
            [~, i1] = min(abs(res_1st{tc}.ta - t_now));
            [~, i2] = min(abs(res_2nd{tc}.ta - t_now));
            ci = clr_overlay(jj,:);

            % 1st order – solid line
            plot(xc, res_1st{tc}.ha(i1,:), '-', ...
                 'LineWidth', lw, 'Color', ci);
            % 2nd order – dashed + sparse marker
            plot(xc, res_2nd{tc}.ha(i2,:), '--', ...
                 'LineWidth', lw, 'Color', ci, ...
                 'Marker', mrk_overlay{jj}, 'MarkerIndices', 1:mrk_every:N, ...
                 'MarkerSize', 5, 'MarkerFaceColor', ci);

            leg_h{end+1} = sprintf('1st, $t=%.2f$ s', t_now);
            leg_h{end+1} = sprintf('2nd, $t=%.2f$ s', t_now);
        end
        hold off; grid on;
        title('$h(x,t)$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        ylabel('$h$',  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        legend(leg_h, 'Interpreter',interp_tex, 'Location','best', ...
               'FontName',fnt, 'FontSize',fs-8, 'NumColumns', 2);
        set(gca, 'FontName',fnt, 'FontSize', fs-4, 'TickLabelInterpreter','latex');

        % ---- subplot 2: q(x,t) ----
        subplot(2,1,2); hold on;
        plot(xc, q0_ic, '-', 'LineWidth', 1.2, 'Color', ic_color);   % IC ghost
        leg_q = {'IC ($t=0$)'};

        for jj = 1:length(t_overlay)
            t_now = t_overlay(jj);
            [~, i1] = min(abs(res_1st{tc}.ta - t_now));
            [~, i2] = min(abs(res_2nd{tc}.ta - t_now));
            ci = clr_overlay(jj,:);

            plot(xc, res_1st{tc}.qa(i1,:), '-', ...
                 'LineWidth', lw, 'Color', ci);
            plot(xc, res_2nd{tc}.qa(i2,:), '--', ...
                 'LineWidth', lw, 'Color', ci, ...
                 'Marker', mrk_overlay{jj}, 'MarkerIndices', 1:mrk_every:N, ...
                 'MarkerSize', 5, 'MarkerFaceColor', ci);

            leg_q{end+1} = sprintf('1st, $t=%.2f$ s', t_now);
            leg_q{end+1} = sprintf('2nd, $t=%.2f$ s', t_now);
        end
        hold off; grid on;
        title('$q(x,t)$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        ylabel('$q=hu$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        legend(leg_q, 'Interpreter',interp_tex, 'Location','best', ...
               'FontName',fnt, 'FontSize',fs-8, 'NumColumns', 2);
        set(gca, 'FontName',fnt, 'FontSize', fs-4, 'TickLabelInterpreter','latex');

        sgtitle(sprintf(['Case (%s): %s --- $\\Delta x = %.4f$, ' ...
                'CFL $= %.2f$, $\\Delta t_{\\mathrm{avg}} \\approx %.4f$ ms'], ...
                cases(tc).tag, cases(tc).name, Delta_x, CFL, dt_avg*1000), ...
                'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs+1);

        exportgraphics(gcf, fullfile('Plots', ...
            sprintf('SWE_overlay_case_%s_dx%.4f.pdf', cases(tc).tag, Delta_x)), ...
            'ContentType','vector');
    end
end

%%  HEATMAP PLOTS
% Two figures: one for 1st order, one for 2nd order.
% Each figure is 2×3: row 1 = h(x,t) for cases a,b,c ; row 2 = q(x,t).
if HeatmapPlots
    order_labels = {'1st order', '2nd order (LW+MC)'};
    order_res    = {res_1st, res_2nd};
    order_tags   = {'1st', '2nd'};

    for order = 1:2
        fig = figure('Units','centimeters','Position',[2 2 42 22], ...
                     'Color','w','PaperPositionMode','auto');

        for tc = 1:3
            R  = order_res{order}{tc};
            Nt = length(R.ta);
            idx = round(linspace(1, Nt, min(Nt,500)));

            % --- row 1: h ---
            subplot(2, 3, tc); 
            pcolor(xc, R.ta(idx), R.ha(idx,:));
            shading interp; colorbar; colormap(gca, parula);
            title(sprintf('(%s) %s --- $h$', cases(tc).tag, cases(tc).name), ...
                  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-4);
            xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-4);
            ylabel('$t$',  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-4);
            set(gca, 'FontName',fnt, 'FontSize', fs-6, 'TickLabelInterpreter','latex');

            % --- row 2: q ---
            subplot(2, 3, 3 + tc);
            pcolor(xc, R.ta(idx), R.qa(idx,:));
            shading interp; colorbar; colormap(gca, parula);
            title(sprintf('(%s) %s --- $q$', cases(tc).tag, cases(tc).name), ...
                  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-4);
            xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-4);
            ylabel('$t$',  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-4);
            set(gca, 'FontName',fnt, 'FontSize', fs-6, 'TickLabelInterpreter','latex');
        end

        dt_avg = mean(diff(order_res{order}{1}.ta));

        sgtitle(sprintf(['%s --- heatmaps, $\\Delta x = %.4f$, ' ...
                'CFL $= %.2f$, $\\Delta t_{\\mathrm{avg}} \\approx %.4f$ ms'], ...
                order_labels{order}, Delta_x, CFL, dt_avg*1000), ...
                'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);

        exportgraphics(gcf, fullfile('Plots', ...
            sprintf('SWE_heatmap_%s_dx%.4f.pdf', order_tags{order}, Delta_x)), ...
            'ContentType','vector');
    end
end

%% CONVERGENCE -- 1st and 2nd order SIDE BY SIDE
if ConvergencePlots
    Ns_conv  = [50, 100, 200, 400, 800];
    t_conv   = 0.5;
    
    % Harmonized color palette (no black)
    clr_conv = [ ...
        0.64 0.08 0.18;   % dark red
        0.85 0.33 0.10;   % orange
        0.00 0.45 0.74;   % blue
        0.47 0.67 0.19;   % green
        0.49 0.18 0.56];  % purple
    
    % Legend entries with N and dx
    leg_entries = arrayfun(@(n) sprintf('$N=%d$, $\\Delta x=%.4f$', n, L/n), ...
                           Ns_conv, 'UniformOutput', false);
    
    % --- h convergence: 1st (left) vs 2nd (right) ---
    fig = figure('Units','centimeters','Position',[2 2 34 16], ...
                 'Color','w','PaperPositionMode','auto');
    
    subplot(1,2,1); hold on;
    subplot(1,2,2); hold on;
    
    for jj = 1:length(Ns_conv)
        Nj  = Ns_conv(jj);
        dxj = L / Nj;
        xj  = ((1:Nj) - 0.5) * dxj;
        h0j = zeros(1,Nj); q0j = zeros(1,Nj);
        h0j(xj < 0.5) = 1.0; h0j(xj >= 0.5) = 0.1;
    
        % 1st order
        [~,~,~, ha1j,~, ta1j] = Godunov1st(h0j, q0j, xj, dxj, g, CFL, t_conv+0.001, [t_conv]);
        [~, i1j] = min(abs(ta1j - t_conv));
        subplot(1,2,1);
        plot(xj, ha1j(i1j,:), 'LineWidth', lw, 'Color', clr_conv(jj,:));
    
        % 2nd order
        [~,~,~, ha2j,~, ta2j] = Godunov2nd(h0j, q0j, xj, dxj, g, CFL, t_conv+0.001, [t_conv]);
        [~, i2j] = min(abs(ta2j - t_conv));
        subplot(1,2,2);
        plot(xj, ha2j(i2j,:), 'LineWidth', lw, 'Color', clr_conv(jj,:));
    end
    
    subplot(1,2,1); hold off; grid on;
    title('1st order --- $h(x,t)$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-2);
    xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    ylabel('$h$',  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    legend(leg_entries, 'Interpreter',interp_tex, 'Location','best', 'FontName',fnt, 'FontSize',fs-8);
    set(gca, 'FontName',fnt, 'FontSize', fs-4, 'TickLabelInterpreter','latex');
    
    subplot(1,2,2); hold off; grid on;
    title('2nd order (LW+MC) --- $h(x,t)$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-2);
    xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    ylabel('$h$',  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    legend(leg_entries, 'Interpreter',interp_tex, 'Location','best', 'FontName',fnt, 'FontSize',fs-8);
    set(gca, 'FontName',fnt, 'FontSize', fs-4, 'TickLabelInterpreter','latex');
    
    % Compute dt_avg from the main run for the sgtitle
    dt_avg_conv = mean(diff(res_1st{1}.ta));

    sgtitle(sprintf(['Case (a): Dam Break --- grid convergence $h$ at $t=%.2f$ s, ' ...
            'CFL$=%.2f$, $\\Delta t_{\\mathrm{avg}} \\approx %.4f$ ms'], ...
            t_conv, CFL, dt_avg_conv*1000), 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    
    exportgraphics(gcf, fullfile('Plots', ...
        sprintf('SWE_convergence_h_caseA_t%.2f.pdf', t_conv)), 'ContentType','vector');
    
    % --- q convergence: 1st (left) vs 2nd (right) ---
    fig = figure('Units','centimeters','Position',[2 2 34 16], ...
                 'Color','w','PaperPositionMode','auto');
    
    subplot(1,2,1); hold on;
    subplot(1,2,2); hold on;
    
    for jj = 1:length(Ns_conv)
        Nj  = Ns_conv(jj);
        dxj = L / Nj;
        xj  = ((1:Nj) - 0.5) * dxj;
        h0j = zeros(1,Nj); q0j = zeros(1,Nj);
        h0j(xj < 0.5) = 1.0; h0j(xj >= 0.5) = 0.1;
    
        % 1st order
        [~,~,~,~, qa1j, ta1j] = Godunov1st(h0j, q0j, xj, dxj, g, CFL, t_conv+0.001, [t_conv]);
        [~, i1j] = min(abs(ta1j - t_conv));
        subplot(1,2,1);
        plot(xj, qa1j(i1j,:), 'LineWidth', lw, 'Color', clr_conv(jj,:));
    
        % 2nd order
        [~,~,~,~, qa2j, ta2j] = Godunov2nd(h0j, q0j, xj, dxj, g, CFL, t_conv+0.001, [t_conv]);
        [~, i2j] = min(abs(ta2j - t_conv));
        subplot(1,2,2);
        plot(xj, qa2j(i2j,:), 'LineWidth', lw, 'Color', clr_conv(jj,:));
    end
    
    subplot(1,2,1); hold off; grid on;
    title('1st order --- $q(x,t)$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-2);
    xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    ylabel('$q=hu$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    legend(leg_entries, 'Interpreter',interp_tex, 'Location','best', 'FontName',fnt, 'FontSize',fs-8);
    set(gca, 'FontName',fnt, 'FontSize', fs-4, 'TickLabelInterpreter','latex');
    
    subplot(1,2,2); hold off; grid on;
    title('2nd order (LW+MC) --- $q(x,t)$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-2);
    xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    ylabel('$q=hu$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    legend(leg_entries, 'Interpreter',interp_tex, 'Location','best', 'FontName',fnt, 'FontSize',fs-8);
    set(gca, 'FontName',fnt, 'FontSize', fs-4, 'TickLabelInterpreter','latex');
    
    sgtitle(sprintf(['Case (a): Dam Break --- grid convergence $q$ at $t=%.2f$ s, ' ...
            'CFL$=%.2f$, $\\Delta t_{\\mathrm{avg}} \\approx %.4f$ ms'], ...
            t_conv, CFL, dt_avg_conv*1000), 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    
    exportgraphics(gcf, fullfile('Plots', ...
        sprintf('SWE_convergence_q_caseA_t%.2f.pdf', t_conv)), 'ContentType','vector');
end


fprintf('\n=== All plots saved to Plots/ folder ===\n');

%% VIDEO: Wave evolution over time
if VideoSol
    vid_fps   = 30;
    vid_dur   = 5;                   % seconds
    n_frames  = vid_fps * vid_dur;   % 150 frames
    
    % Info string for video titles
    info_str = sprintf('$\\Delta x = %.4f$, CFL $= %.2f$, $N = %d$, $\\Delta t_{\\mathrm{avg}} \\approx %.4f$ ms', ...
                       Delta_x, CFL, N, mean(diff(res_1st{1}.ta))*1000);
    
    for tc = 1:3
        R1 = res_1st{tc};
        R2 = res_2nd{tc};
    
        % Pick n_frames evenly spaced time indices from the full history
        Nt1 = length(R1.ta);
        Nt2 = length(R2.ta);
        idx1 = round(linspace(1, Nt1, n_frames));
        idx2 = round(linspace(1, Nt2, n_frames));
    
        % Global y-axis limits (fixed across all frames for smooth animation)
        h_min_plot = min([min(R1.ha(:)), min(R2.ha(:))]) - 0.02;
        h_max_plot = max([max(R1.ha(:)), max(R2.ha(:))]) + 0.02;
        q_min_plot = min([min(R1.qa(:)), min(R2.qa(:))]) - 0.02;
        q_max_plot = max([max(R1.qa(:)), max(R2.qa(:))]) + 0.02;
    
        % Create video file
        vid_name = fullfile('Plots', sprintf('SWE_video_case_%s.mp4', cases(tc).tag));
        v = VideoWriter(vid_name, 'MPEG-4');
        v.FrameRate = vid_fps;
        v.Quality   = 95;
        open(v);
    
        fig = figure('Units','pixels','Position',[100 100 1280 720], ...
                     'Color','w','Visible','off');
    
        for k = 1:n_frames
            t1 = R1.ta(idx1(k));
    
            % ---- h(x,t) ----
            subplot(2,1,1); cla; hold on;
            plot(xc, R1.h0, '-', 'LineWidth', 1.0, 'Color', ic_color);   % IC ghost
            plot(xc, R1.ha(idx1(k),:), '-',  'LineWidth', lw, 'Color', [0.0 0.45 0.74]);
            plot(xc, R2.ha(idx2(k),:), '--', 'LineWidth', lw, 'Color', [0.85 0.33 0.10]);
            hold off; grid on;
            ylim([h_min_plot, h_max_plot]);
            xlim([0, L]);
            title(sprintf('Case (%s): %s --- $h(x,t)$ at $t = %.3f$ s \\quad [%s]', ...
                  cases(tc).tag, cases(tc).name, t1, info_str), ...
                  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-2);
            xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
            ylabel('$h$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
            legend({'IC ($t=0$)','1st order','2nd order (LW+MC)'}, ...
                   'Interpreter',interp_tex, 'Location','northeast', ...
                   'FontName',fnt, 'FontSize',fs-6);
            set(gca, 'FontName',fnt, 'FontSize',fs-4, 'TickLabelInterpreter','latex');
    
            % ---- q(x,t) ----
            subplot(2,1,2); cla; hold on;
            plot(xc, R1.q0, '-', 'LineWidth', 1.0, 'Color', ic_color);   % IC ghost
            plot(xc, R1.qa(idx1(k),:), '-',  'LineWidth', lw, 'Color', [0.0 0.45 0.74]);
            plot(xc, R2.qa(idx2(k),:), '--', 'LineWidth', lw, 'Color', [0.85 0.33 0.10]);
            hold off; grid on;
            ylim([q_min_plot, q_max_plot]);
            xlim([0, L]);
            title(sprintf('$q(x,t)$ at $t = %.3f$ s', t1), ...
                  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-2);
            xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
            ylabel('$q = hu$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
            legend({'IC ($t=0$)','1st order','2nd order (LW+MC)'}, ...
                   'Interpreter',interp_tex, 'Location','northeast', ...
                   'FontName',fnt, 'FontSize',fs-6);
            set(gca, 'FontName',fnt, 'FontSize',fs-4, 'TickLabelInterpreter','latex');
    
            % Capture frame
            drawnow;
            frame = getframe(fig);
            writeVideo(v, frame);
        end
    
        close(v);
        close(fig);
        fprintf('Video saved: %s\n', vid_name);
    end
    fprintf('\n=== All videos saved to Plots/ folder ===\n');
end