function [h_snaps, q_snaps, t_snaps_out, h_all, q_all, t_all] = ...
    Godunov1stLog(h0, q0, xc, Delta_x, g, CFL, Tf, snap_times, h_ref)
% -------------------------------------------------------------------------
% First-order Godunov scheme for the LOG-PRESSURE SWE.
%
% Differences from Godunov1st.m:
%   1. CFL uses c = sqrt(g*(1 + log(h/h_ref))) instead of sqrt(g*h)
%   2. Calls RiemannSolverLog (with h_ref) instead of RiemannSolver
%
% INPUT/OUTPUT: same as Godunov1st, plus h_ref.
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
        % c(h) = sqrt(g*(1 + log(h/h_ref))) instead of sqrt(g*h)
        u_vel = q ./ h;
        w     = 1 + log(h / h_ref);
        w     = max(w, 1e-10);          % safety for hyperbolicity
        c     = sqrt(g * w);
        max_s = max(abs(u_vel) + c);
        dt    = CFL * Delta_x / max_s;

        if t + dt > Tf, dt = Tf - t; end

        % Ghost cells (reflecting wall BC — unchanged)
        h_ext = [h(1), h, h(N)];
        q_ext = [-q(1), q, -q(N)];

        flux_h = zeros(1, N+1);
        flux_q = zeros(1, N+1);
        for i = 1:N+1
            hL = h_ext(i);   hR = h_ext(i+1);
            uL = q_ext(i)/hL; uR = q_ext(i+1)/hR;

            % ---- LOG-PRESSURE Riemann solver ----
            [fh, fq] = RiemannSolverLog(hL, uL, hR, uR, g, h_ref);

            flux_h(i) = fh; flux_q(i) = fq;
        end

        h = h - (dt/Delta_x) * (flux_h(2:end) - flux_h(1:end-1));
        q = q - (dt/Delta_x) * (flux_q(2:end) - flux_q(1:end-1));
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
