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

% PLOT OF NUMERICAL SOLUTION
figure;

% ----- SUBPLOT 1: magnitude -----
subplot(2,1,1);
plot(x_h, abs(uh_full), '-o', 'LineWidth', 1.5);
grid on;
ylabel('|u_h(x)|');
title(sprintf('Magnitude of u_h (h = %g,  nEl = %g)', femregion.h, nEl_calc));

% ----- SUBPLOT 2: phase -----
subplot(2,1,2);
plot(x_h, angle(uh_full), '-o', 'LineWidth', 1.5);
grid on;
ylabel('\angle u_h(x)');
xlabel('x');
title('Phase of u_h');
title(sprintf('Phase of u_h (h = %g,  nEl = %g)', femregion.h, nEl_calc));

% ----- SAVE IMAGE -----
% if Data.save_sol_images
%     filename = sprintf('solution3_nEL_%g.png', nEl_calc);
%     exportgraphics(gcf, fullfile('Plots', filename), 'Resolution', 300);
% end


% PLOT OF NUMERICAL SOLUTION VS EXACT
x = linspace(Data.domain(1),Data.domain(2),1000);
uex_full = Data.uex(x, Data.omega, Data.ro, Data.vel);

figure;

fs = 25;   % SAME font size as your other figures

% ----- SUBPLOT 1: magnitude -----
subplot(2,1,1);
plot(x_h, abs(uh_full), '-o', 'LineWidth', 1.5); hold on;
plot(x, abs(uex_full), '--', 'LineWidth', 1.8);
grid on;

ylabel('$|u_h(x)|$', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman', ...
       'FontSize',fs);

xlabel('$x$', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman', ...
       'FontSize',fs);

legend({'$u_h$ (numerical)', '$u$ (exact)'}, ...
       'Interpreter','latex', ...
       'Location','best', ...
       'FontName','Times New Roman', ...
       'FontSize',fs);

title(sprintf('Magnitude of $u_h$ ($h=%.3g$, $n_{\\mathrm{El}}=%d$)', ...
      femregion.h, nEl_calc), ...
      'Interpreter','latex', ...
      'FontName','Times New Roman', ...
      'FontSize',fs);

% ----- SUBPLOT 2: phase -----
subplot(2,1,2);
plot(x_h, angle(uh_full), '-o', 'LineWidth', 1.8); hold on;
plot(x, angle(uex_full), '--', 'LineWidth', 1.8);
grid on;

ylabel('$\angle u_h(x)$', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman', ...
       'FontSize',fs);

xlabel('$x$', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman', ...
       'FontSize',fs);

legend({'$u_h$ (numerical)', '$u$ (exact)'}, ...
       'Interpreter','latex', ...
       'Location','best', ...
       'FontName','Times New Roman', ...
       'FontSize',fs);

title(sprintf('Phase of $u_h$ ($h=%.3g$, $n_{\\mathrm{El}}=%d$)', ...
      femregion.h, nEl_calc), ...
      'Interpreter','latex', ...
      'FontName','Times New Roman', ...
      'FontSize',fs);

% ----- SAVE IMAGE -----
% if Data.save_sol_images
%     filename = sprintf('solution3_nEL_withExact_%g.png', nEl_calc);
%     exportgraphics(gcf, fullfile('Plots', filename), 'Resolution', 300);
% end
if Data.save_sol_images
    filename = sprintf('solution5_nEL_withExact_%g.pdf', nEl_calc);
    exportgraphics(gcf, fullfile('Plots', filename), ...
                   'ContentType','vector');
end

end
