function [Lh, Lq] = MUSCLLog(h, q, N, Delta_x, g, h_min, h_ref)
% -------------------------------------------------------------------------
% MUSCL spatial operator for the LOG-PRESSURE SWE.
% Identical to MUSCL.m except it calls RiemannSolverLog with h_ref.
%
% INPUT:
%   - h, q          Cell averages
%   - N             Number of cells
%   - Delta_x       Cell size
%   - g             Gravity
%   - h_min         Positivity tolerance
%   - h_ref         Reference depth for log-pressure
% OUTPUT:
%   - Lh, Lq        MUSCL operator components
% -------------------------------------------------------------------------

    % BC: extended with 2 ghost cells per side (reflecting wall)
    h_e = [h(2), h(1), h, h(N), h(N-1)];
    q_e = [-q(2), -q(1), q, -q(N), -q(N-1)];

    % MC-limited slopes (unchanged — limiter is flux-independent)
    sh = zeros(1, N);
    sq = zeros(1, N);
    for i = 1:N
        ii = i + 2;
        sh(i) = MCslope(h_e(ii-1), h_e(ii), h_e(ii+1), Delta_x);
        sq(i) = MCslope(q_e(ii-1), q_e(ii), q_e(ii+1), Delta_x);
    end

    % Lax-Wendroff reconstruction at interfaces
    hL_int = zeros(1, N+1);  qL_int = zeros(1, N+1);
    hR_int = zeros(1, N+1);  qR_int = zeros(1, N+1);

    for i = 1:N+1
        if i == 1
            hL_int(i) = h(1);           qL_int(i) = -q(1);
            hR_int(i) = h(1) - 0.5*Delta_x*sh(1);
            qR_int(i) = q(1) - 0.5*Delta_x*sq(1);
        elseif i == N+1
            hL_int(i) = h(N) + 0.5*Delta_x*sh(N);
            qL_int(i) = q(N) + 0.5*Delta_x*sq(N);
            hR_int(i) = h(N);           qR_int(i) = -q(N);
        else
            hL_int(i) = h(i-1) + 0.5*Delta_x*sh(i-1);
            qL_int(i) = q(i-1) + 0.5*Delta_x*sq(i-1);
            hR_int(i) = h(i)   - 0.5*Delta_x*sh(i);
            qR_int(i) = q(i)   - 0.5*Delta_x*sq(i);
        end
    end

    hL_int = max(hL_int, h_min);
    hR_int = max(hR_int, h_min);

    % Riemann solves with LOG-PRESSURE flux
    flux_h = zeros(1, N+1);
    flux_q = zeros(1, N+1);
    for i = 1:N+1
        uL = qL_int(i) / hL_int(i);
        uR = qR_int(i) / hR_int(i);
        [fh, fq] = RiemannSolverLog(hL_int(i), uL, hR_int(i), uR, g, h_ref);
        flux_h(i) = fh;  flux_q(i) = fq;
    end

    Lh = -(1/Delta_x) * (flux_h(2:end) - flux_h(1:end-1));
    Lq = -(1/Delta_x) * (flux_q(2:end) - flux_q(1:end-1));
end
