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

% if Data.snapshot
% 
%     %======================================================================
%     % VIDEO SETUP
%     %======================================================================
% 
%     outFolder = 'Plots';
%     if ~exist(outFolder, 'dir')
%         mkdir(outFolder);
%     end
% 
%     videoFile = fullfile(outFolder, [Data.name, '.mp4']);
% 
%     v = VideoWriter(videoFile, 'MPEG-4');
%     v.FrameRate = 30;          % playback speed
%     open(v);
% 
%     %======================================================================
%     % SNAPSHOT / MOVIE LOOP
%     %======================================================================
% 
%     for k = 1:Data.frameSkip:size(eta_snap,2)
% 
%         t = (k-1) * Data.dt;
% 
%         Snapshot(Femregion, eta_snap(:,k), q_snap(:,k), Data, t);
% 
%         drawnow;
% 
%         % Capture and write frame
%         frame = getframe(gcf);
%         writeVideo(v, frame);
% 
%     end
% 
%     %======================================================================
%     % FINALIZE VIDEO
%     %======================================================================
% 
%     close(v);
% end

% Solution at certaint time
Solutions = PostProcessing(Data, Femregion, eta_snap(:,round(Data.PicTime/Data.dt)), q_snap(:,round(Data.PicTime/Data.dt)));

% Surface plots
if Data.surf
    time = 0:Data.dt:Data.T;

    % x coords consistent with reduced solution length
    x = Femregion.coord(:,1);
    if length(x) ~= size(eta_snap,1)
        x = x(1:size(eta_snap,1));
    end

    [tt,xx] = meshgrid(time, x);

    figure(3); clf;
    
    fs = 25;   % SAME font size as all other figures
    
    surf(tt, xx, eta_snap, 'EdgeColor', 'none');
    grid on;
    
    xlabel('$t$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);
    
    ylabel('$x$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);
    
    zlabel('$\eta(x,t)$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);

    title(sprintf('Surface plot of $\\eta$  with $\\Delta t=%.2f$\\,ms, $h=%.2f$', ...
          Data.PicTime, Data.dt*1000, 1/Femregion.ne), ...
          'Interpreter','latex', ...
          'FontName','Times New Roman', ...
          'FontSize',fs);
    
    
    if Data.save_sol_images
        filename = sprintf('SurfacePlot_ETA_%s_nEL_%g_p_%g_dt_%g.pdf', ...
                           Data.name, Femregion.ne, Data.p, Data.dt);
        exportgraphics(gcf, fullfile('Plots', filename), ...
                       'ContentType','vector');
    end


    figure(4); clf;

    surf(tt, xx, q_snap, 'EdgeColor', 'none');
    grid on;
    
    xlabel('$t$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);
    
    ylabel('$x$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);
    
    zlabel('$q(x,t)$', ...
           'Interpreter','latex', ...
           'FontName','Times New Roman', ...
           'FontSize',fs);
    
    title(sprintf('Surface plot of $q$  with $\\Delta t=%.2f$\\,ms, $h=%.2f$', ...
          Data.PicTime, Data.dt*1000, 1/Femregion.ne), ...
          'Interpreter','latex', ...
          'FontName','Times New Roman', ...
          'FontSize',fs);
    
    if Data.save_sol_images
        filename = sprintf('SurfacePlot_q_%s_nEL_%g_p_%g_dt_%g.pdf', ...
                           Data.name, Femregion.ne, Data.p, Data.dt);
        exportgraphics(gcf, fullfile('Plots', filename), ...
                       'ContentType','vector');
    end

end

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
Errors = [];
if (Data.calc_errors)
    [Errors] = ComputeErrors(Data,Femregion,Solutions);
end



