function [errors] = ComputeErrors(Data, femregion, solutions)

% ---- periodic: expand reduced vectors back to full length for connectivity
if strcmp(Data.boundary,'PP') && length(solutions.eta_h) == femregion.ndof-1
    eta_h_full = [solutions.eta_h; solutions.eta_h(1)];
    q_h_full   = [solutions.q_h;   solutions.q_h(1)];
else
    eta_h_full = solutions.eta_h;
    q_h_full   = solutions.q_h;
end

% =========================
% ETA errors
% =========================
Data_eta = Data;
Data_eta.uex     = Data.etaex;
Data_eta.graduex = Data.gradetaex;

[E_eta_L2, E_eta_SEMI_H1] = Calc_L2_H1_errors(femregion, eta_h_full, Data_eta);

% For inf-norm use reduced vectors directly (matching what you plot)
err_eta = solutions.eta_h - solutions.eta_ex;

% =========================
% Q errors
% =========================
Data_q = Data;
Data_q.uex     = Data.qex;
Data_q.graduex = Data.gradqex;

[E_q_L2, E_q_SEMI_H1] = Calc_L2_H1_errors(femregion, q_h_full, Data_q);

err_q = solutions.q_h - solutions.q_ex;

% =========================
% pack output
% =========================
errors = struct( ...
    'eta_L2',   E_eta_L2, ...
    'eta_H1_S', E_eta_SEMI_H1, ...
    'eta_H1',   sqrt(E_eta_L2.^2 + E_eta_SEMI_H1.^2), ...
    'eta_inf',  norm(full(err_eta), inf), ...
    'q_L2',     E_q_L2, ...
    'q_H1_S',   E_q_SEMI_H1, ...
    'q_H1',     sqrt(E_q_L2.^2 + E_q_SEMI_H1.^2), ...
    'q_inf',    norm(full(err_q), inf) );
end
