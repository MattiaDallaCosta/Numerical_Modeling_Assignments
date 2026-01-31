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

figure;
plot(x_h, uh_full, '-o', 'LineWidth', 1.8); hold on;

% Evaluate exact solution at the DOFs
uex_full = Data.uex(x, Data.omega, Data.ro, Data.vel);
plot(x, uex_full, '--', 'LineWidth', 1.8);

grid on;
legend('u_h (numerical)', 'u (exact)', 'Location', 'best', fontsize=18);

title(sprintf('u_h solution for h = %g and nEl = %g', femregion.h, nEl_calc), fontsize=17);
xlabel('x', fontsize=17);
ylabel('u(x)', fontsize=17);

if (Data.save_sol_images)
    filename = sprintf('solution2_nEL_%g.png', nEl_calc);
    exportgraphics(gcf, fullfile('Plots', filename), 'Resolution', 300);
end
end
