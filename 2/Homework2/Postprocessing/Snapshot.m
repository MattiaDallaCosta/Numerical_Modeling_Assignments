function Snapshot(femregion, eta_h, q_h, Data, t)
%% Snapshot(femregion, eta_h, q_h, Data, t)
%==========================================================================
% Plot eta_h and q_h at time t (optionally with exact solutions).
%==========================================================================

x1 = femregion.domain(1,1);
x2 = femregion.domain(1,2);

% Coordinates (must match the reduced system if periodic)
x = femregion.coord(:,1);

% If your periodic reduction removed last dof, then eta_h/q_h are shorter:
% just plot with x(1:end-1)
if length(eta_h) ~= length(x)
    x = x(1:length(eta_h));
end

% Exact solutions if provided
has_etaex = isfield(Data,'etaex') && isa(Data.etaex,'function_handle');
has_qex   = isfield(Data,'qex')   && isa(Data.qex,'function_handle');

eta_ex = [];
q_ex   = [];
if has_etaex, eta_ex = Data.etaex(x,t); end
if has_qex,   q_ex   = Data.qex(x,t);   end

% Axis limits (auto, but keep your style if you prefer fixed)
M1 = max(full(eta_h)); m1 = min(full(eta_h));
M2 = max(full(q_h));   m2 = min(full(q_h));
pad1 = 0.05*(M1-m1 + eps); pad2 = 0.05*(M2-m2 + eps);

figure(1); clf;

subplot(2,1,1);
plot(x, full(eta_h), 'LineWidth', 1.5); hold on;
if has_etaex
    plot(x, eta_ex, '--', 'LineWidth', 1.5);
    legend('\eta_h','\eta_{ex}','Location','best');
else
    legend('\eta_h','Location','best');
end
title(['\eta(x,t) at time t = ', num2str(t), ' s']);
xlabel('x'); ylabel('\eta');
axis([x1, x2, m1-pad1, M1+pad1]);

subplot(2,1,2);
plot(x, full(q_h), 'LineWidth', 1.5); hold on;
if has_qex
    plot(x, q_ex, '--', 'LineWidth', 1.5);
    legend('q_h','q_{ex}','Location','best');
else
    legend('q_h','Location','best');
end
title(['q(x,t) at time t = ', num2str(t), ' s']);
xlabel('x'); ylabel('q');
axis([x1, x2, m2-pad2, M2+pad2]);

drawnow;
pause(0.000000000002);
end
