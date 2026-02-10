function [uh] = PlotSol(femregion, Data, uh)
%% Snapshot(femregion, uh, string)
%==========================================================================
% PLOT THE NUMERICAL SOLUTION u_h ON THE DOFS
%==========================================================================
%    called in MainHMZ.m
%
%    INPUT:
%          femregion   : (struct)  see CreateFemregion.m
%          uh          : (sparse(ndof,1) real/complex) solution vector
%          string      : (string, opzionale) extra text for the title
%

% convert sparse vector to full
uh_full = full(uh);

x_h = femregion.coord(:,1);
x = linspace(Data.domain(1),Data.domain(2),1000);
nEl_calc = (Data.domain(2) - Data.domain(1))/femregion.h;

fs        = 30;
lw        = 1.9;
fnt       = 'Times New Roman';
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

fig1 = figure('Units','centimeters','Position',[2 2 30 16], ...
              'Color','w','PaperPositionMode','auto');
hold on; box on; grid on;
set(gca, 'GridLineStyle','-', 'GridAlpha',0.15);

plot(x_h, uh_full, '-o', 'LineWidth',lw, 'MarkerSize',6, 'Color',clr(1,:));

% Evaluate exact solution at the DOFs
uex_full = Data.uex(x, Data.omega, Data.ro, Data.vel);
plot(x, uex_full, '--', 'LineWidth',lw, 'Color',clr(2,:));

legend({'$u_h$ (numerical)', '$u$ (exact)'}, ...
       'Interpreter',interp_tex, ...
       'Location','best', ...
       'FontName',fnt, ...
       'FontSize',fs-8, ...
       'Box','on');

title(sprintf('$u_h$ solution for $h = %g$ and $n_{\\mathrm{El}} = %g$', ...
      femregion.h, nEl_calc), ...
      'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
xlabel('$x$',    'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
ylabel('$u(x)$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);

set(gca, 'FontName',fnt, 'FontSize',fs-4, 'TickLabelInterpreter','latex');

if Data.save_sol_images
    filename = sprintf('solution2_nEL_%d.pdf', nEl_calc);
    exportgraphics(fig1, fullfile('Plots', filename), 'ContentType','vector');
end


% figure;
% plot(x_h, uh_full, '-o', 'LineWidth', 1.8); hold on;
% 
% % Evaluate exact solution at the DOFs
% uex_full = Data.uex(x, Data.omega, Data.ro, Data.vel);
% plot(x, uex_full, '--', 'LineWidth', 1.8);
% 
% grid on;
% 
% fs = 25;   % tuned for 0.45\linewidth inclusion
% 
% legend({'$u_h$ (numerical)', '$u$ (exact)'}, ...
%        'Interpreter','latex', ...
%        'Location','best', ...
%        'FontName','Times New Roman', ...
%        'FontSize',fs);
% 
% title(sprintf('u_h solution for h = %g and nEl = %g', ...
%       femregion.h, nEl_calc), ...
%       'FontName','Times New Roman', ...
%       'FontSize',fs);
% 
% xlabel('x', 'FontName','Times New Roman', 'FontSize',fs);
% ylabel('u(x)', 'FontName','Times New Roman', 'FontSize',fs);
% 
% % if (Data.save_sol_images)
% %     filename = sprintf('solution2_nEL_%g.png', nEl_calc);
% %     exportgraphics(gcf, fullfile('Plots', filename), 'Resolution', 300);
% % end
% if Data.save_sol_images
%     filename = sprintf('solution2_nEL_%d.pdf', nEl_calc);
%     exportgraphics(gcf, fullfile('Plots', filename), ...
%                    'ContentType','vector');
% end
end
