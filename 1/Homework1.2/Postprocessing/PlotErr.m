function PlotErr(hVec, eVecL2, eVecH1, pL2, pH1, Data)
%% PlotErr(hVec, eVecL2, eVecH1, pL2, pH1)
%==========================================================================
% Plot L2 and H1-seminorm errors vs h in log-log scale
% with reference slopes h^{pL2} and h^{pH1}
%
% INPUT:
%   hVec   : vector of mesh sizes h (e.g. [fem1.h, fem2.h, fem3.h, fem4.h])
%   eVecL2 : vector of L2 errors   (e.g. [err1.L2, err2.L2, err3.L2, err4.L2])
%   eVecH1 : vector of H1 errors   (e.g. [err1.H1, err2.H1, err3.H1, err4.H1])
%   pL2    : (optional) expected L2 convergence rate (default = 2)
%   pH1    : (optional) expected H1 convergence rate (default = 1)
%==========================================================================

if nargin < 4 || isempty(pL2)
    pL2 = 2;   % expected L2 order
end
if nargin < 5 || isempty(pH1)
    pH1 = 1;   % expected H1 order
end

% Build reference lines scaled to match last error point
C_L2 = eVecL2(end) / (hVec(end)^pL2);
C_H1 = eVecH1(end) / (hVec(end)^pH1);

refL2 = C_L2 * hVec.^pL2;
refH1 = C_H1 * hVec.^pH1;

figure;

%---------------------- L2 error ----------------------%
%hs1 = subplot(2,1,1);
loglog(hVec, refL2, '-+b', 'LineWidth', 2); hold on;
loglog(hVec, eVecL2, '-or', 'LineWidth', 2);
grid on;

% optional: fit experimental slope
pfitL2 = polyfit(log(hVec), log(eVecL2), 1);
rateL2 = pfitL2(1);

legend( ...
    sprintf('C h^{%d}', pL2), ...
    sprintf('||u - u_h||_{L^2}, rate ≈ %.2f', rateL2), ...
    'Location', 'best', fontsize=18);

ylabel('L^2-error', fontsize=17);
xlabel('h', fontsize=17);
title('L^2 error vs h', fontsize=18);
hs1.FontSize = 12;

%---------------------- H1 error ----------------------%
% hs2 = subplot(2,1,2);
% loglog(hVec, refH1, '-+b', 'LineWidth', 2); hold on;
% loglog(hVec, eVecH1, '-or', 'LineWidth', 2);
% grid on;
% 
% pfitH1 = polyfit(log(hVec), log(eVecH1), 1);
% rateH1 = pfitH1(1);
% 
% legend( ...
%     sprintf('C h^{%d}', pH1), ...
%     sprintf('||u - u_h||_{H^1}, rate ≈ %.2f', rateH1), ...
%     'Location', 'best');
% 
% ylabel('H^1-error (seminorm)');
% xlabel('h');
% title('H^1 error vs h');
% hs2.FontSize = 12;

exportgraphics(gcf, fullfile('Plots', 'Error.png'), 'Resolution', 300);

end
