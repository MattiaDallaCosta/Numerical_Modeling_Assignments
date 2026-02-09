function [h_snaps, q_snaps, t_snaps_out, h_all, q_all, t_all] = Godunov1st(h0, q0, xc, Delta_x, g, CFL, Tf, snap_times)
% -------------------------------------------------------------------------
% First-order Godunov finite volume scheme for the 1D SWE.
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
    h_min = 1e-6; % tolerance for min h (h must be kept positive)

    
    % Create buffers with max time sample 200000
    max_steps = 200000;
    h_buf = zeros(max_steps, N);    % 200000 X N (time X space)
    q_buf = zeros(max_steps, N);    % 200000 X N (time X space)
    t_buf = zeros(max_steps, 1);    % 200000 X 1 (time)

    % Save initial conditions in the buffer
    step = 1;

    h_buf(step,:) = h;
    q_buf(step,:) = q;
    t_buf(step) = 0;

    % Prepare buffers to save time snapshot
    n_snaps = length(snap_times);       % Number of snaps
    h_snaps = zeros(n_snaps, N);        % nsnaps X N
    q_snaps = zeros(n_snaps, N);        % nsnaps X N
    t_snaps_out = zeros(1, n_snaps);    % 1 X nsnaps
    snap_capt = false(1, n_snaps);      % Buffer to save if snaps are saved

    % Capture t = 0 snapshot if requested
    for k = 1:n_snaps
        if snap_times(k) == 0
            h_snaps(k,:) = h;
            q_snaps(k,:) = q;
            t_snaps_out(k) = 0;
            snap_capt(k) = true;
        end
    end

    % Starting time
    t = 0;

    % Iterate t = t + Deltat until we reach final time. Deltat is defined
    % thanks to CFL condition
    while t < Tf - 1e-14

        % Enforce positivity of h
        h = max(h, h_min);

        % Compute velocities:
        u_vel = q ./ h; 
        c     = sqrt(g * h);
        max_s = max(abs(u_vel) + c);

        % Comute dt using CFL condition
        dt = CFL * Delta_x / max_s;
        % make sure we don't exceed final time
        if t + dt > Tf
            dt = Tf - t;
        end

        % Ghost cells for reflecting wall BC:
        %   Left wall:  h_ghost = h(1),   q_ghost = -q(1)
        %   Right wall: h_ghost = h(N),   q_ghost = -q(N)
        % The continuity of h across the boundary ensures no spatial
        % variation of depth across the interface
        h_ext = [h(1), h, h(N)];
        % The sign flip on q mirrors the moemntum, enforcing u=0 at walls.
        q_ext = [-q(1), q, -q(N)];

        % Preallocate to store  Godunov numerical flux at each interface 
        % i+1/2, i=0,...,N
        flux_h = zeros(1, N+1);
        flux_q = zeros(1, N+1);

        for i = 1:N+1
            % Compute, for each interface, h_L, h_R, q_L and q_R.
            hL = h_ext(i);
            hR = h_ext(i+1);
            uL = q_ext(i)/hL;
            uR = q_ext(i+1)/hR;

            % Solve the Riemann problem 
            [fh, fq] = RiemannSolver(hL, uL, hR, uR, g);

            % Store solution
            flux_h(i) = fh;
            flux_q(i) = fq;
        end

        % Use the FV scheme to find h at new time step
        h = h - (dt/Delta_x) * (flux_h(2:end) - flux_h(1:end-1));
        q = q - (dt/Delta_x) * (flux_q(2:end) - flux_q(1:end-1));

        % Update time
        t = t + dt;

        % Save solutions at the of time step
        step = step + 1;
        if step <= max_steps
            h_buf(step,:) = h;  % h at time t
            q_buf(step,:) = q;  % q at time t
            t_buf(step) = t;    % time t
        end

        % Capture requested snapshots
        for k = 1:n_snaps
            if ~snap_capt(k) && t >= snap_times(k) - 1e-10
                h_snaps(k,:) = h;
                q_snaps(k,:) = q;
                t_snaps_out(k) = t;
                snap_capt(k) = true;
            end
        end
    end

    % h, q and t buffers contain all solutions. Let's store at maximum
    % max_steps = 200000 time steps.
    actual = min(step, max_steps);
    h_all = h_buf(1:actual,:);
    q_all = q_buf(1:actual,:);
    t_all = t_buf(1:actual);
end
