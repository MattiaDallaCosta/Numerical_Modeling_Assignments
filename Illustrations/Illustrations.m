%%  Plot Linear (hat) basis functions and Lagrange polynomials
clear; close all; clc;

% Formatting (consistent with homework style)
fs   = 22;
lw   = 1.9;
fnt  = 'Times New Roman';
interp_tex = 'latex';

% Color palette (standard, solid, non-pastel)
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

Nh = 10;                        % number of elements
L  = 1;                         % domain length
h  = L / Nh;                    % element size
Npts = Nh + 1;                  % number of nodes (including endpoints)
x_nodes = linspace(0, L, Npts); % node positions
x_fine  = linspace(0, L, 1000); % fine grid for plotting

fig1 = figure('Units','centimeters','Position',[2 2 30 16], ...
              'Color','w','PaperPositionMode','auto');
hold on; box on; grid on;
set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.15);

leg_entries = cell(1, Npts);

for j = 1:Npts
    phi = zeros(size(x_fine));
    xj = x_nodes(j);

    for k = 1:length(x_fine)
        xx = x_fine(k);
        if j > 1 && xx >= x_nodes(j-1) && xx <= xj
            % Left leg: element k_{j-1} = [x_{j-1}, x_j]
            phi(k) = (xx - x_nodes(j-1)) / h;
        elseif j < Npts && xx >= xj && xx <= x_nodes(j+1)
            % Right leg: element k_j = [x_j, x_{j+1}]
            phi(k) = (x_nodes(j+1) - xx) / h;
        elseif xx == xj
            phi(k) = 1;
        end
    end

    ci = clr(mod(j-1, size(clr,1)) + 1, :);
    plot(x_fine, phi, '-', 'LineWidth', lw, 'Color', ci);
    leg_entries{j} = sprintf('$\\varphi_{%d}(x)$', j-1);
end

hold off;
xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
ylabel('$\varphi_j(x)$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
xlim([0 L]); ylim([-0.05 1.2]);
legend(leg_entries, 'Interpreter',interp_tex, 'FontName',fnt, ...
       'FontSize',fs-8, 'Location','eastoutside', 'Box','on', 'NumColumns',1);
set(gca, 'FontName',fnt, 'FontSize',fs-4, 'TickLabelInterpreter','latex');

sgtitle(sprintf('Linear basis functions, $N_h = %d$, $L = %d$, $N_{\\mathrm{pts}} = %d$', ...
        Nh, L, Npts), 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);

exportgraphics(fig1, 'linear_basis_functions.pdf', 'ContentType','vector');
fprintf('Saved: linear_basis_functions.pdf\n');

% Lagrange polynomials on GLL points, degree p = 1..4

% GLL nodes for each polynomial degree (precomputed, well-known values)
gll_nodes = { ...
    [-1, 1], ...                                     % p = 1
    [-1, 0, 1], ...                                  % p = 2
    [-1, -sqrt(1/5), sqrt(1/5), 1], ...              % p = 3
    [-1, -sqrt(3/7), 0, sqrt(3/7), 1]};             % p = 4

xi_fine = linspace(-1, 1, 500);

% Colors for each degree
clr_deg = [ ...
    0.00  0.45  0.74;   % blue
    0.85  0.33  0.10;   % orange
    0.47  0.67  0.19;   % green
    0.64  0.08  0.18];  % dark red

% Line styles for basis functions within each subplot
ls_basis = {'-', '--', '-.', ':', '-'};
ls_basis = {'-', '-', '-', '-', '-'};

fig2 = figure('Units','centimeters','Position',[2 2 34 24], ...
              'Color','w','PaperPositionMode','auto');

for p = 1:4
    subplot(2, 2, p);
    hold on; box on; grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.25);

    nodes = gll_nodes{p};
    n_nodes = length(nodes);
    leg_p = cell(1, n_nodes + 1);

    % Plot each Lagrange basis polynomial l_j(xi)
    for j = 1:n_nodes
        lj = ones(size(xi_fine));
        for m = 1:n_nodes
            if m ~= j
                lj = lj .* (xi_fine - nodes(m)) / (nodes(j) - nodes(m));
            end
        end

        ci = clr(mod(j-1, size(clr,1)) + 1, :);
        style = ls_basis{mod(j-1, length(ls_basis)) + 1};
        plot(xi_fine, lj, style, 'LineWidth', lw, 'Color', ci);
        leg_p{j} = sprintf('$\\varphi_{%d}(\\xi)$', j);
    end

    % Mark GLL nodes on x-axis
    plot(nodes, zeros(size(nodes)), 'k^', 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.3 0.3 0.3], 'LineWidth', 1.2);
    leg_p{end} = 'GLL nodes';

    % Zero line
    plot([-1 1], [0 0], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);

    hold off;
    xlabel('$\xi$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-4);
    ylabel('$\varphi_j(\xi)$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-4);
    title(sprintf('Degree $p = %d$ \\quad ($%d$ GLL nodes)', p, n_nodes), ...
          'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-4);
    xlim([-1.05 1.05]);
    ylim([-0.4 1.15]);
    legend(leg_p, 'Interpreter',interp_tex, 'FontName',fnt, ...
           'FontSize',fs-10, 'Location','best', 'Box','on');
    set(gca, 'FontName',fnt, 'FontSize',fs-6, 'TickLabelInterpreter','latex');
end

sgtitle('Lagrange polynomials on Gauss--Lobatto--Legendre nodes, $\xi \in [-1,\,1]$', ...
        'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);

exportgraphics(fig2, 'gll_lagrange_polynomials.pdf', 'ContentType','vector');
fprintf('Saved: gll_lagrange_polynomials.pdf\n');
%% Finite Volume Method Discretization 

% Parameters
N = 10;                         % number of cells
L = 1;                          % domain length
dx = L / N;                     % cell width
x_faces = linspace(0, L, N+1);  % cell interfaces x_{i-1/2}, x_{i+1/2}
x_centers = 0.5*(x_faces(1:end-1) + x_faces(2:end)); % cell centers

% Continuous function to approximate
f = @(x) 0.15 + 0.6*exp(-12*(x - 0.35).^2) + 0.25*sin(2.5*pi*x);
x_fine = linspace(0, L, 500);
y_fine = f(x_fine);

% Integral's cell averages
u_avg = zeros(1, N);
for i = 1:N
    xq = linspace(x_faces(i), x_faces(i+1), 200);
    u_avg(i) = trapz(xq, f(xq)) / dx;
end

% Figure
fig = figure();
set(fig, 'Color', 'w');
ax = axes('Position', [0.08 0.13 0.88 0.82]);
hold on; box on; grid off;

% Color
col_func  = [0.0 0.0 0.0];        
col_step  = [0.0 0.45 0.74];     
col_fill  = [0.0 0.45 0.74 0.12]; 
col_face  = [0.5 0.5 0.5];       
col_center = [0.85 0.33 0.10];  

% Draw cell averages as filled rectangles
for i = 1:N
    % Filled rectangle
    patch([x_faces(i) x_faces(i+1) x_faces(i+1) x_faces(i)], ...
          [0 0 u_avg(i) u_avg(i)], ...
          col_step, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    % Top horizontal line of the step
    plot([x_faces(i) x_faces(i+1)], [u_avg(i) u_avg(i)], ...
         '-', 'Color', col_step, 'LineWidth', 1.8);
    % Vertical lines at interfaces (step jumps)
    if i > 1
        plot([x_faces(i) x_faces(i)], [u_avg(i-1) u_avg(i)], ...
             '-', 'Color', col_step, 'LineWidth', 1.8);
    else
        plot([x_faces(i) x_faces(i)], [0 u_avg(i)], ...
             '-', 'Color', col_step, 'LineWidth', 1.8);
    end
end

% Last vertical line
plot([x_faces(end) x_faces(end)], [0 u_avg(end)], ...
     '-', 'Color', col_step, 'LineWidth', 1.8);

% Highlight one cell
i_highlight = 5;
x_L = x_faces(i_highlight);
x_R = x_faces(i_highlight+1);
h_val = u_avg(i_highlight);

% Plot continuous function
plot(x_fine, y_fine, '-', 'Color', col_func, 'LineWidth', 2.0);

% Cell interfaces (dashed vertical lines)
for i = 1:N+1
    plot([x_faces(i) x_faces(i)], [0 max(y_fine)*1.08], ...
         '--', 'Color', col_face, 'LineWidth', 0.6);
end

% Grid points at cell centers
plot(x_centers, zeros(size(x_centers)), 'o', ...
     'MarkerSize', 5, 'MarkerFaceColor', col_center, ...
     'MarkerEdgeColor', col_center, 'LineWidth', 1.0);

% Cell interface markers on x-axis
plot(x_faces, zeros(size(x_faces)), '|', ...
     'MarkerSize', 8, 'Color', col_face, 'LineWidth', 1.2);

% Annotations

% Label highlighted cell
mid_x = x_centers(i_highlight);
text(mid_x, h_val + 0.04, ['$\bar{u}_{', num2str(i_highlight), '}$'], ...
     'Interpreter', 'latex', 'FontSize', 11, 'Color', col_step, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Label adjacent cells
text(x_centers(i_highlight-1), u_avg(i_highlight-1) + 0.04, ...
     ['$\bar{u}_{', num2str(i_highlight-1), '}$'], ...
     'Interpreter', 'latex', 'FontSize', 10, 'Color', col_step, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(x_centers(i_highlight+1), u_avg(i_highlight+1) + 0.04, ...
     ['$\bar{u}_{', num2str(i_highlight+1), '}$'], ...
     'Interpreter', 'latex', 'FontSize', 10, 'Color', col_step, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Label cell centers on x-axis
for idx = [i_highlight-1, i_highlight, i_highlight+1]
    text(x_centers(idx), -0.045, ['$x_{', num2str(idx), '}$'], ...
         'Interpreter', 'latex', 'FontSize', 9, 'Color', col_center, ...
         'HorizontalAlignment', 'center');
end

% Label cell interfaces
text(x_faces(i_highlight), -0.075, ['$x_{', num2str(i_highlight), '-\frac{1}{2}}$'], ...
     'Interpreter', 'latex', 'FontSize', 9, 'Color', col_face, ...
     'HorizontalAlignment', 'center');
text(x_faces(i_highlight+1), -0.075, ['$x_{', num2str(i_highlight), '+\frac{1}{2}}$'], ...
     'Interpreter', 'latex', 'FontSize', 9, 'Color', col_face, ...
     'HorizontalAlignment', 'center');

% Delta x brace annotation
y_brace = -0.11;
plot([x_faces(i_highlight) x_faces(i_highlight+1)], [y_brace y_brace], ...
     '-', 'Color', col_step, 'LineWidth', 1.2);
plot([x_faces(i_highlight) x_faces(i_highlight)], [y_brace-0.01 y_brace+0.01], ...
     '-', 'Color', col_step, 'LineWidth', 1.2);
plot([x_faces(i_highlight+1) x_faces(i_highlight+1)], [y_brace-0.01 y_brace+0.01], ...
     '-', 'Color', col_step, 'LineWidth', 1.2);
text(mid_x, y_brace - 0.03, '$\Delta x$', ...
     'Interpreter', 'latex', 'FontSize', 10, 'Color', col_step, ...
     'HorizontalAlignment', 'center');

% Legend
h1 = plot(NaN, NaN, '-', 'Color', col_func, 'LineWidth', 2.0);
h2 = plot(NaN, NaN, '-', 'Color', col_step, 'LineWidth', 1.8);
h3 = plot(NaN, NaN, 'o', 'MarkerSize', 5, 'MarkerFaceColor', col_center, ...
     'MarkerEdgeColor', col_center);
legend([h1 h2 h3], {'$u(x)$', '$\bar{u}_i$ (cell average)', 'Grid points'}, ...
       'Interpreter', 'latex', 'FontSize', 15, 'Location', 'best', ...
       'Box', 'on');

% Axes
set(ax, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 12);
xlim([-0.02 L+0.02]);
ylim([-0.16 max(y_fine)*1.15]);


set(ax, 'Position', [0.08 0.13 0.88 0.78]);  % shrink axes slightly
tit = sgtitle('Finite volume approximation', 'Interpreter','latex', 'FontSize',16);

% Export
exportgraphics(fig, 'fvm_discretization.pdf', 'ContentType', 'vector');
fprintf('Saved: fvm_discretization.pdf\n');

%% Comparison: Standard vs Logarithmic Flux

clear; close all; clc;

% Parameters
g     = 9.81;
h_ref = 0.05;
h     = linspace(h_ref/exp(1) + 1e-4, 1.2, 500); % avoid c_log becoming imaginary

% Standard SWE
c_std = sqrt(g * h);
P_std = 0.5 * g * h.^2;

% Logarithmic flux
w     = 1 + log(h / h_ref);
c_log = sqrt(g * w);
P_log = g * h .* log(h / h_ref);

% Colors
col_std = [0.0  0.45 0.74];   
col_log = [0.85 0.33 0.10];   
col_ref = [0.5  0.5  0.5 ];   

% Figure
fig = figure('Units','centimeters','Position',[2 2 18 8]);
set(fig, 'Color', 'w');

% ─── Subplot 1: Wave speed ───
ax1 = subplot(1,2,1);
hold on; box on; grid on;
set(ax1, 'GridLineStyle', ':', 'GridAlpha', 0.3);

plot(h, c_std, '-',  'Color', col_std, 'LineWidth', 1.8);
plot(h, c_log, '--', 'Color', col_log, 'LineWidth', 1.8);


xlabel('$h$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$c(h)$', 'Interpreter', 'latex', 'FontSize', 13);
title('Wave speed', 'Interpreter', 'latex', 'FontSize', 13);
legend({'$c = \sqrt{g\,h}$', ...
        '$c = \sqrt{g\left(1+\log \frac{h}{h_{\mathrm{ref}}}\right)}$'}, ...
       'Interpreter', 'latex', 'FontSize', 9, 'Location', 'northwest', 'Box', 'on');
set(ax1, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
xlim([0 1.2]);
ylim([0 max(c_log)*1.08]);

% ─── Subplot 2: Pressure P(h) ───
ax2 = subplot(1,2,2);
hold on; box on; grid on;
set(ax2, 'GridLineStyle', ':', 'GridAlpha', 0.3);

plot(h, P_std, '-',  'Color', col_std, 'LineWidth', 1.8);
plot(h, P_log, '--', 'Color', col_log, 'LineWidth', 1.8);


xlabel('$h$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$P(h)$', 'Interpreter', 'latex', 'FontSize', 13);
title('Pressure law', 'Interpreter', 'latex', 'FontSize', 13);
legend({'$P = \frac{1}{2}\,g\,h^2$', ...
        '$P = g\,h\,\log\!\left(\frac{h}{h_{\mathrm{ref}}}\right)$'}, ...
       'Interpreter', 'latex', 'FontSize', 13, 'Location', 'northwest', 'Box', 'on');
set(ax2, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
xlim([0 1.2]);
ylim([min(P_log)*1.1 max(P_log)*1.08]);

% Title
sgtitle(['Wave speed and pressure: standard vs logarithmic flux ($h_{\mathrm{ref}} = ', ...
         num2str(h_ref), '$)'], ...
        'Interpreter', 'latex', 'FontSize', 15);

% Export
exportgraphics(fig, 'flux_comparison.pdf', 'ContentType', 'vector', 'Resolution', 300);
fprintf('Saved: flux_comparison.pdf\n');
