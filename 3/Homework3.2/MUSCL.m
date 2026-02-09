
function [Lh, Lq] = MUSCL(h, q, N, Delta_x, g, h_min)
% -------------------------------------------------------------------------
% Computes the MUSCL (Monotone Upstream-centred Schemes for Conservation
% Laws) spatial discretization operator for the 1D SWE.
%
% Given cell averages (h, q), this function:
%   1. Applies reflecting wall ghost cells
%   2. Computes MC-limited slopes for each component
%   3. Builds reconstructed interface states (U^-, U^+)
%   4. Solves exact Riemann problems at each interface
%   5. Returns L(U) = -(1/Delta_x) * (F_{i+1/2} - F_{i-1/2})
%
% INPUT:
%   - h             depth
%   - q             momentum
%   - N             spatial mesh
%   - Delta_x       spatial step
%   - g             gravity
%   - h_min           positivity of h
% OUTPUT:
%   - Lh            MUSCL operator first component
%   - Lq            MUSCL operator second component
% -------------------------------------------------------------------------

    % BC. Bigger estention of boundaries to apply L-W scheme
    h_e = [h(2), h(1), h, h(N), h(N-1)];
    q_e = [-q(2), -q(1), q, -q(N), -q(N-1)];

    % MC-limited slopes for h and q
    sh = zeros(1, N);
    sq = zeros(1, N);
    for i = 1:N
        ii = i + 2;  % offset for ghost cells
        sh(i) = MCslope(h_e(ii-1), h_e(ii), h_e(ii+1), Delta_x);
        sq(i) = MCslope(q_e(ii-1), q_e(ii), q_e(ii+1), Delta_x);
    end

    % LAX-WELDROFF second order scheme
    %   U^-_{i+1/2} = U_i     + (1/2)*Delta_x * sigma_i
    %   U^+_{i+1/2} = U_{i+1} - (1/2)*Delta_x * sigma_{i+1}
    hL_int = zeros(1, N+1);
    qL_int = zeros(1, N+1);
    hR_int = zeros(1, N+1);
    qR_int = zeros(1, N+1);

    for i = 1:N+1
        if i == 1
            % Left boundary (reflecting wall)
            hL_int(i) = h(1);
            qL_int(i) = -q(1);
            hR_int(i) = h(1) - 0.5*Delta_x*sh(1);
            qR_int(i) = q(1) - 0.5*Delta_x*sq(1);
        elseif i == N+1
            % Right boundary (reflecting wall)
            hL_int(i) = h(N) + 0.5*Delta_x*sh(N);
            qL_int(i) = q(N) + 0.5*Delta_x*sq(N);
            hR_int(i) = h(N);
            qR_int(i) = -q(N);
        else
            % Interior interface
            hL_int(i) = h(i-1) + 0.5*Delta_x*sh(i-1);
            qL_int(i) = q(i-1) + 0.5*Delta_x*sq(i-1);
            hR_int(i) = h(i)   - 0.5*Delta_x*sh(i);
            qR_int(i) = q(i)   - 0.5*Delta_x*sq(i);
        end
    end

    % Enforce positivity
    hL_int = max(hL_int, h_min);
    hR_int = max(hR_int, h_min);

    % Solve exact Riemann problem at each interface (as Godunov1)
    flux_h = zeros(1, N+1);
    flux_q = zeros(1, N+1);
    for i = 1:N+1
        uL = qL_int(i) / hL_int(i);
        uR = qR_int(i) / hR_int(i);
        [fh, fq] = RiemannSolver(hL_int(i), uL, hR_int(i), uR, g);
        flux_h(i) = fh;  flux_q(i) = fq;
    end

    % Return L(U) = -(1/Delta_x) * (F_{i+1/2} - F_{i-1/2})
    Lh = -(1/Delta_x) * (flux_h(2:end) - flux_h(1:end-1));
    Lq = -(1/Delta_x) * (flux_q(2:end) - flux_q(1:end-1));
end

