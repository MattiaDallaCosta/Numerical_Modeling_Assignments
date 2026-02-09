%% Finite Volume Method Discretization 
clear; close all; clc;

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
