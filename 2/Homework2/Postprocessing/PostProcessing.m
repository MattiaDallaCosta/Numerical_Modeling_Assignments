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
    eta_ex_T = Data.etaex(x, Data.PicTime);
end
if isfield(Data,'qex') && isa(Data.qex,'function_handle')
    q_ex_T = Data.qex(x, Data.PicTime);
end

% Plot final time
if Data.visual_graph
    figure(2); clf;
    
    fs = 25;   % SAME font size as all other figures
    
    % ---------- SUBPLOT 1: eta ----------
    subplot(2,1,1);
    plot(x, full(eta_h_T), 'LineWidth', 1.9); hold on;
    if ~isempty(eta_ex_T)
        plot(x, eta_ex_T, '--', 'LineWidth', 1.9);
    end
    grid on;
    
    title(sprintf('$\\eta$ at time %.1f\\,s with $\\Delta t=%.2f$\\,ms, $h=%.2f$', ...
          Data.PicTime, Data.dt*1000, 1/femregion.ne), ...
          'Interpreter','latex', ...
          'FontName','Times New Roman', ...
          'FontSize',fs);
    
    xlabel('$x$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);
    
    ylabel('$\eta$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);
    
    if ~isempty(eta_ex_T)
        legend({'$\eta_h$', '$\eta_{\mathrm{ex}}$'}, ...
               'Interpreter','latex', ...
               'Location','best', ...
               'FontName','Times New Roman', ...
               'FontSize',fs);
    end
    
    
    % ---------- SUBPLOT 2: q ----------
    subplot(2,1,2);
    plot(x, full(q_h_T), 'LineWidth', 1.9); hold on;
    if ~isempty(q_ex_T)
        plot(x, q_ex_T, '--', 'LineWidth', 1.9);
    end
    grid on;
    
    title(sprintf('$q$ at time %.1f\\,s with $\\Delta t=%.2f$\\,ms, $h=%.2f$', ...
          Data.PicTime, Data.dt*1000, 1/femregion.ne), ...
          'Interpreter','latex', ...
          'FontName','Times New Roman', ...
          'FontSize',fs);
    
    xlabel('$x$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);
    
    ylabel('$q$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);
    
    if ~isempty(q_ex_T)
        legend({'$q_h$', '$q_{\mathrm{ex}}$'}, ...
               'Interpreter','latex', ...
               'Location','best', ...
               'FontName','Times New Roman', ...
               'FontSize',fs);
    end


    if Data.save_sol_images
        filename = sprintf('%s_nEL_%g_p_%g_dt_%g_t_%g.pdf', Data.name, femregion.ne, Data.p, Data.dt, Data.PicTime);
        exportgraphics(gcf, fullfile('Plots', filename), ...
                   'ContentType','vector');
    end
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
