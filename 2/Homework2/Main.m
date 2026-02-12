function [Errors,Solutions,Femregion,Data] = Main(Data, nEl)
%%
%    INPUT:
%          Data    : (struct) Data struct
%          nEl     : (int)    Number of mesh elements
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) finite element space
%
%          Data        : (struct)  Data struct
%

fprintf('============================================================\n')
fprintf(['Solving test ', Data.name, ' with ',num2str(nEl),' elements \n']);

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = CreateMesh(Data,nEl);

% Inside this function, for each sub-element of the domain, the xwlgl
% function is called. Such function returns the weigths and the nodes of 
% the corresponding Legendre Gauss-Lobatto quadrature formula in the  
% interval [a,b]

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[Femregion] = CreateFemregion(Data,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES
%==========================================================================

[D_no_bc,M_no_bc] = Matrix1D(Data,Femregion);

% In Matrix1d the basis_and_der_at_lgl function is used to obtain the value 
% of the basis functions and of  the spectral Legendre Gauss Lobatto 
% derivative matrix d at the np LGL nodes x (on [-1,1]). Thanks to that the
% M and D matrix are calculated:
%   D_ij = ∫ phi_j * (phi_i)_x dx
%   M_ij = ∫ phi_i * phi_j dx

% Observe that if the polynomial degree is p, and we have nEl number of
% elements in the domain. The total basis functions are:
%       nEl * (p + 1) - (nEl - 1) <<--- (subtract the common nodes)
% Indeed the matrix are of dimension (nEl * p + 1, nEl * p + 1)

% Apply Boundary Conditions
D = D_no_bc;
M = M_no_bc;

if strcmp(Data.boundary,'PP')

    % Periodic boundary conditions reduce dof(1) by one. We indeed add the
    % last row and column to the first, in order to diminish the matric of
    % 1 dimension

    % Merge end DOF into first DOF (rows and columns)
    D(1,:) = D(1,:) + D(end,:);
    D(:,1) = D(:,1) + D(:,end);

    M(1,:) = M(1,:) + M(end,:);
    M(:,1) = M(:,1) + M(:,end);

    % Remove the redundant last DOF
    D = D(1:end-1, 1:end-1);
    M = M(1:end-1, 1:end-1);
end

%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================

[b_nbc] = Rhs1D(Data,Femregion);
b = b_nbc;

% BC
if strcmp(Data.boundary,'PP')
    b(1,:) = b(1,:) + b(end,:);
    b = b(1:end-1,:);
end

% ============================================================
% BUILD THE FULL BLOCK MATRICES
% ============================================================
N  = size(M,1);
Z  = sparse(N,N);

Mfull = [M  Z;
         Z  M];

Sfull = [Z      -D;
        -Data.g*Data.H*D   Z];

% Forcing is usually placed only in first equation, we indeed set second
% block to zero
Bfull = [b; sparse(N, size(b,2))];

% ============================================================
% FRICTION HANDLER
% ============================================================

% Add friction: q_t + gH eta_x = -gamma q => add +gamma*M*q in the second 
% equation

if isfield(Data,'gamma') && Data.gamma > 0
    Sfull = Sfull + [Z, Z;
                     Z, Data.gamma * M];
end

%==========================================================================
% BUILD INITIAL CONDITIONS AND TIME INTEGRATION (theta-method)
%==========================================================================

% Theta definition
theta = 1/2;

% Time definition
time = 0:Data.dt:Data.T;
nT   = numel(time);

% Initial condition vectors
x = Femregion.coord(:,1);

eta0 = Data.eta0(x);
q0   = Data.q0(x);

% Periodic reduction for ICs
if strcmp(Data.boundary,'PP')
    eta0 = eta0(1:end-1);
    q0   = q0(1:end-1);
end

% Build matrix to store U = [ eta q]^T for each time sample
U = zeros(2*N, nT);
U(:,1) = [eta0; q0];

% Let's rewrite the theta-method bringing all k+1 steps on the left and k
% steps on the right:
%   (M/dt + theta S) U_{k+1} = (M/dt - (1-theta) S) U_k
%         A    *     U_{k+1} =        G       *     U_k
A = (1/Data.dt) * Mfull + theta * Sfull;
G = (1/Data.dt) * Mfull - (1-theta) * Sfull;

% Reflecting wall BC (q=0 at x=0,L): prepare modified system matrix
A_bc = A;
if strcmp(Data.boundary,'RW')
    idx_q_left  = N + 1; % Points of the first node
    idx_q_right = 2*N;   % Points of the last node

    A_bc(idx_q_left, :) = 0;
    A_bc(idx_q_left, idx_q_left) = 1;

    A_bc(idx_q_right, :) = 0;
    A_bc(idx_q_right, idx_q_right) = 1;
end

% For reach time step find U in k+1 by solving A \ rhs of the equation
for k = 1:nT-1
    rhs = G * U(:,k) + theta * Bfull(:,k+1) + (1-theta) * Bfull(:,k);
    %rhs = G * U(:,k);
    if strcmp(Data.boundary,'RW')
        rhs(idx_q_left)  = 0;
        rhs(idx_q_right) = 0;
        U(:,k+1) = A_bc \ rhs;
    else
        U(:,k+1) = A \ rhs;
    end
end

% Split solution again
eta_snap = U(1:N, :);
q_snap   = U(N+1:end, :);

%==========================================================================
% POST-PROCESSING
%==========================================================================

% Movie / snapshots
if Data.snapshot
    t = 0;
    % for i = 1:size(eta_snap,2)
    %     Snapshot(Femregion, eta_snap(:,i), q_snap(:,i), Data, t);
    %     t = t + Data.dt;
    % end
    for k = 1:Data.frameSkip:size(eta_snap,2)
        t = (k-1)*Data.dt;
        Snapshot(Femregion, eta_snap(:,k), q_snap(:,k), Data, t);
    end
end

%======================================================================
% VIDEO: Wave evolution over time (eta & q)
%======================================================================
x_ex = linspace(Data.domain(1),Data.domain(2),1000);
if Data.snapshot

    % ---------------------------
    % Video parameters
    % ---------------------------
    vid_fps  = 30;

    % If user provided a duration, use it; otherwise infer from T
    if isfield(Data,'vid_dur') && ~isempty(Data.vid_dur)
        vid_dur = Data.vid_dur;                 % seconds
    else
        vid_dur = min(10, max(2, Data.T));      % fallback (2–10 s)
    end

    n_frames = max(1, round(vid_fps * vid_dur));

    % Frame picking from available time history (use Data.frameSkip if present)
    Nt = size(eta_snap, 2);

    % If you want to respect Data.frameSkip first, build an "available list"
    if isfield(Data,'frameSkip') && Data.frameSkip > 0
        avail = 1:Data.frameSkip:Nt;
    else
        avail = 1:Nt;
    end

    % Now pick n_frames evenly spaced indices from avail
    idx = round(linspace(1, numel(avail), min(n_frames, numel(avail))));
    idx = avail(idx);

    % Avoid duplicates just in case
    idx = unique(idx, 'stable');
    n_frames = numel(idx);

    % ---------------------------
    % Axes limits (fixed)
    % ---------------------------
    eta_min_plot = min(eta_snap(:)) - 0.02;
    eta_max_plot = max(eta_snap(:)) + 0.02;
    q_min_plot   = min(q_snap(:))   - 0.02;
    q_max_plot   = max(q_snap(:))   + 0.02;

    % ---------------------------
    % Output folder and filename
    % ---------------------------
    outFolder = 'Plots';
    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    % Add discretization info in filename (optional)
    vid_name = fullfile(outFolder, sprintf('%s_nEl_%d_p_%d_dt_%g.mp4', ...
                    Data.name, Femregion.ne, Data.p, Data.dt));

    v = VideoWriter(vid_name, 'MPEG-4');
    v.FrameRate = vid_fps;
    v.Quality   = 95;
    open(v);

    % ---------------------------
    % Styling (reuse your globals below if you prefer)
    % ---------------------------
    fs         = 28;
    lw         = 2.0;
    fnt        = 'Times New Roman';
    interp_tex = 'latex';

    x = Femregion.coord(:,1);
    if strcmp(Data.boundary,'PP')
        % periodic reduction has size N already, x must match snapshots
        if length(x) ~= size(eta_snap,1)
            x = x(1:size(eta_snap,1));
        end
    else
        if length(x) ~= size(eta_snap,1)
            x = x(1:size(eta_snap,1));
        end
    end

    % Domain limits
    xL = min(x);
    xR = max(x);

    % Info string (title suffix)
    info_str = sprintf('$nEl=%d$, $p=%d$, $\\Delta t=%.3f$ ms, $T=%.2f$ s', ...
                       Femregion.ne, Data.p, Data.dt*1000, Data.T);

    fig = figure('Units','pixels','Position',[100 100 1280 720], ...
                 'Color','w','Visible','off');

    % ---------------------------
    % Frame loop
    % ---------------------------
    for kk = 1:n_frames
        k = idx(kk);
        t = (k-1) * Data.dt;

        % ---- eta(x,t) ----
        subplot(2,1,1); cla; hold on;
        plot(x, eta_snap(:,k), '-', 'LineWidth', lw);
        plot(x_ex, Data.etaex(x_ex,t), '--', 'LineWidth', lw);
        legend('$\eta_h$','$\eta_{ex}$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        %legend('$\eta_h$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        hold off; grid on;
        xlim([xL, xR]);
        ylim([eta_min_plot, eta_max_plot]);
        title(sprintf('Linearized SWE. $\\eta(x,t)$ at $t=%.3f$ s \\quad [%s]', t, info_str), ...
              'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-2);
        xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        ylabel('$\eta$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        set(gca, 'FontName',fnt, 'FontSize',fs-6, 'TickLabelInterpreter','latex');

        % ---- q(x,t) ----
        subplot(2,1,2); cla; hold on;
        plot(x, q_snap(:,k), '-', 'LineWidth', lw);
        plot(x_ex, Data.qex(x_ex,t), '--', 'LineWidth', lw);
        hold off; grid on;
        xlim([xL, xR]);
        ylim([q_min_plot, q_max_plot]);
        title(sprintf('$q(x,t)$ at $t=%.3f$ s', t), ...
              'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs-2);
        legend('$q_h$','$q_{ex}$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        %legend('$q_h$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        xlabel('$x$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        ylabel('$q$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
        set(gca, 'FontName',fnt, 'FontSize',fs-6, 'TickLabelInterpreter','latex');

        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end

    close(v);
    close(fig);
    fprintf('Video saved: %s\n', vid_name);
end

% Solution at certaint time
Solutions = PostProcessing(Data, Femregion, eta_snap(:,round(Data.PicTime/Data.dt)), q_snap(:,round(Data.PicTime/Data.dt)));

% Surface plots
fs         = 35;
lw         = 1.9;
fnt        = 'Times New Roman';
interp_tex = 'latex';

clr = [ ...
    0.00  0.45  0.74;   %  1  blue
    0.85  0.33  0.10;   %  2  orange
    0.47  0.67  0.19;   %  3  green
    0.49  0.18  0.56;   %  4  purple
    0.93  0.69  0.13;   %  5  yellow
    0.30  0.75  0.93;   %  6  cyan
    0.64  0.08  0.18;   %  7  dark red
    0.25  0.25  0.25;   %  8  dark gray
    0.00  0.60  0.50;   %  9  teal
    0.75  0.00  0.75;   % 10  magenta
    0.60  0.60  0.60];  % 11  light gray

if Data.surf

    time = 0:Data.dt:Data.T;
    x    = Femregion.coord(:,1);
    if length(x) ~= size(eta_snap,1)
        x = x(1:size(eta_snap,1));
    end

    % ---- Optional decimation for very large meshes ----
    maxPts = 400;                       % max grid points per axis
    step_x = max(1, floor(length(x)    / maxPts));
    step_t = max(1, floor(length(time) / maxPts));
    x_d    = x(1:step_x:end);
    t_d    = time(1:step_t:end);
    eta_d  = eta_snap(1:step_x:end, 1:step_t:end);
    q_d    = q_snap(1:step_x:end,   1:step_t:end);
    [tt, xx] = meshgrid(t_d, x_d);

    % ---------- Surface: eta ----------
    fig2a = figure('Units','centimeters','Position',[2 2 30 20], ...
                   'Color','w','PaperPositionMode','auto');
    surf(tt, xx, eta_d, 'EdgeColor','none');
    box on; grid on;
    set(gca, 'GridLineStyle','-', 'GridAlpha',0.15);
    colorbar;

    xlabel('$t$',          'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    ylabel('$x$',          'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    zlabel('$\eta(x,t)$',  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    title(sprintf('Surface plot of $\\eta$ with $\\Delta t=%.2f$\\,ms, $h=%.2f$, $p=%g.$', ...
          Data.dt*1000, 1/Femregion.ne, Data.p), ...
          'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    set(gca, 'FontName',fnt, 'FontSize',fs-4, 'TickLabelInterpreter','latex');

    if Data.save_sol_images
        filename = sprintf('SurfacePlot_ETA_%s_nEL_%g_p_%g_dt_%g.pdf', ...
                           Data.name, Femregion.ne, Data.p, Data.dt);
        exportgraphics(fig2a, fullfile('Plots', filename), 'Resolution',300);
    end

    % ---------- Surface: q ----------
    fig2b = figure('Units','centimeters','Position',[2 2 30 20], ...
                   'Color','w','PaperPositionMode','auto');
    surf(tt, xx, q_d, 'EdgeColor','none');
    box on; grid on;
    set(gca, 'GridLineStyle','-', 'GridAlpha',0.15);
    colorbar;

    xlabel('$t$',       'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    ylabel('$x$',       'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    zlabel('$q(x,t)$',  'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    title(sprintf('Surface plot of $q$ with $\\Delta t=%.2f$\\,ms, $h=%.2f$, $p=%g.$', ...
          Data.dt*1000, 1/Femregion.ne, Data.p), ...
          'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    set(gca, 'FontName',fnt, 'FontSize',fs-4, 'TickLabelInterpreter','latex');

    if Data.save_sol_images
        filename = sprintf('SurfacePlot_q_%s_nEL_%g_p_%g_dt_%g.pdf', ...
                           Data.name, Femregion.ne, Data.p, Data.dt);
        exportgraphics(fig2b, fullfile('Plots', filename), 'Resolution',1000);
    end
end

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
Errors = [];
if (Data.calc_errors)
    [Errors] = ComputeErrors(Data,Femregion,Solutions);
end



