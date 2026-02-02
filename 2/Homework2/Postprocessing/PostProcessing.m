function [solutions] = PostProcessing(Data, femregion, eta_h_T, q_h_T)
%% PostProcessing(Data,femregion,eta_h_T,q_h_T)
%==========================================================================
% Post processing at final time T for eta and q.
%==========================================================================

fprintf('\n Plot the solution ... \n');

x = femregion.dof(:,1);

% Handle periodic reduction mismatch
if length(eta_h_T) ~= length(x)
    x = x(1:length(eta_h_T));
end

% Exact solutions at final time if available
eta_ex_T = [];
q_ex_T   = [];
if isfield(Data,'etaex') && isa(Data.etaex,'function_handle')
    eta_ex_T = Data.etaex(x, Data.T);
end
if isfield(Data,'qex') && isa(Data.qex,'function_handle')
    q_ex_T = Data.qex(x, Data.T);
end

% Plot final time
if Data.visual_graph
    figure(2); clf;

    subplot(2,1,1);
    plot(x, full(eta_h_T), 'LineWidth', 1.5); hold on;
    if ~isempty(eta_ex_T), plot(x, eta_ex_T, '--', 'LineWidth', 1.5); end
    title('\eta at final time T'); xlabel('x'); ylabel('\eta');
    if ~isempty(eta_ex_T), legend('\eta_h','\eta_{ex}','Location','best'); end

    subplot(2,1,2);
    plot(x, full(q_h_T), 'LineWidth', 1.5); hold on;
    if ~isempty(q_ex_T), plot(x, q_ex_T, '--', 'LineWidth', 1.5); end
    title('q at final time T'); xlabel('x'); ylabel('q');
    if ~isempty(q_ex_T), legend('q_h','q_{ex}','Location','best'); end
end

% Save
solutions = struct( ...
    'x', x, ...
    'eta_h', eta_h_T, ...
    'q_h', q_h_T, ...
    'eta_ex', eta_ex_T, ...
    'q_ex', q_ex_T );

fprintf('============================================================\n')
end
