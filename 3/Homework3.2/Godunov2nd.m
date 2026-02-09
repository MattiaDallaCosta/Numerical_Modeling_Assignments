function [h_snaps, q_snaps, t_snaps_out, h_all, q_all, t_all] = ...
    Godunov2nd(h0, q0, xc, Delta_x, g, CFL, Tf, snap_times)
% -------------------------------------------------------------------------
% Second-order Godunov finite volume scheme for the 1D SWE.
%
% INPUT:
%   - h0            depth IC
%   - q0            momentum IC
%   - xc            centers of each cell
%   - Delta_x       spatial step
%   - g             gravity
%   - CFL           CFL parameter
%   - Tf            Final time
%   - snap_times    Time at which store solution for plots
% OUTPUT:
%   - h_snaps       Store solution h for certaint times (for snapshots)
%   - q_snaps       Store solution q for certaint times (for snapshots)
%   - t_snaps_out   Store relative time of the snapshots
%   - h_all         Store all h(x,t) solutions for heatmap
%   - q_all         Store all q(x,t) solutions for heatmap
%   - t_all         Store all times for heatmap
% -------------------------------------------------------------------------

    N = length(h0);
    h = h0;
    q = q0;
    h_min = 1e-6;

    max_steps = 200000;
    h_buf = zeros(max_steps, N);
    q_buf = zeros(max_steps, N);
    t_buf = zeros(max_steps, 1);
    h_buf(1,:) = h;
    q_buf(1,:) = q;
    t_buf(1) = 0;
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

        % CFL
        u_vel = q ./ h;
        c     = sqrt(g * h);
        max_s = max(abs(u_vel) + c);
        dt    = CFL * Delta_x / max_s;

        if t + dt > Tf, dt = Tf - t; end

        % RK2 STAGE 1, I want to compute the intermidiate step U* = (h*  q*)^T 
        % with: U* = U^n + dt * L(U^n).
        
        % The MUSCL operator L is evaluated by the following function.
        [Lh_n, Lq_n] = MUSCL(h, q, N, Delta_x, g, h_min);

        % Store U* components
        h_star = h + dt * Lh_n;
        q_star = q + dt * Lq_n;
        h_star = max(h_star, h_min); % ensure h* positivity

        % RK2 STAGE 2, the next time step solution can be computed by the 
        % following scheme: U^{n+1} = 1/2 U^n + 1/2 (U* + dt*L(U*))
        
        % MUSCL operator computation.
        [Lh_s, Lq_s] = MUSCL(h_star, q_star, N, Delta_x, g, h_min);

        % Save h and q
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