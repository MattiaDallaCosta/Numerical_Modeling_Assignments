function [h_snaps, q_snaps, t_snaps_out, h_all, q_all, t_all] = ...
    Godunov2ndLog(h0, q0, xc, Delta_x, g, CFL, Tf, snap_times, h_ref)
% -------------------------------------------------------------------------
% Second-order SSP-RK2 Godunov scheme for the LOG-PRESSURE SWE.
%
% Differences from Godunov2nd.m:
%   1. CFL uses c = sqrt(g*(1 + log(h/h_ref)))
%   2. Calls MUSCLLog (which calls RiemannSolverLog with h_ref)
%
% INPUT/OUTPUT: same as Godunov2nd, plus h_ref.
% -------------------------------------------------------------------------

    N = length(h0);
    h = h0;  q = q0;
    h_min = 1e-6;

    max_steps = 200000;
    h_buf = zeros(max_steps, N);
    q_buf = zeros(max_steps, N);
    t_buf = zeros(max_steps, 1);
    h_buf(1,:) = h; q_buf(1,:) = q; t_buf(1) = 0;
    step = 1;

    n_snaps = length(snap_times);
    h_snaps = zeros(n_snaps, N);
    q_snaps = zeros(n_snaps, N);
    t_snaps_out = zeros(1, n_snaps);
    snap_capt = false(1, n_snaps);
    for k = 1:n_snaps
        if snap_times(k) == 0
            h_snaps(k,:) = h; q_snaps(k,:) = q;
            t_snaps_out(k) = 0; snap_capt(k) = true;
        end
    end

    t = 0;
    while t < Tf - 1e-14

        h = max(h, h_min);

        % ---- CFL with LOG-PRESSURE celerity ----
        u_vel = q ./ h;
        w     = 1 + log(h / h_ref);
        w     = max(w, 1e-10);
        c     = sqrt(g * w);
        max_s = max(abs(u_vel) + c);
        dt    = CFL * Delta_x / max_s;

        if t + dt > Tf, dt = Tf - t; end

        % SSP-RK2 Stage 1:  U* = U^n + dt * L(U^n)
        [Lh_n, Lq_n] = MUSCLLog(h, q, N, Delta_x, g, h_min, h_ref);
        h_star = h + dt * Lh_n;
        q_star = q + dt * Lq_n;
        h_star = max(h_star, h_min);

        % SSP-RK2 Stage 2:  U^{n+1} = 1/2 U^n + 1/2 (U* + dt*L(U*))
        [Lh_s, Lq_s] = MUSCLLog(h_star, q_star, N, Delta_x, g, h_min, h_ref);
        h = 0.5 * h + 0.5 * (h_star + dt * Lh_s);
        q = 0.5 * q + 0.5 * (q_star + dt * Lq_s);
        h = max(h, h_min);

        t = t + dt;

        step = step + 1;
        if step <= max_steps
            h_buf(step,:) = h; q_buf(step,:) = q; t_buf(step) = t;
        end

        for k = 1:n_snaps
            if ~snap_capt(k) && t >= snap_times(k) - 1e-10
                h_snaps(k,:) = h; q_snaps(k,:) = q;
                t_snaps_out(k) = t; snap_capt(k) = true;
            end
        end
    end

    actual = min(step, max_steps);
    h_all = h_buf(1:actual,:);
    q_all = q_buf(1:actual,:);
    t_all = t_buf(1:actual);
end
