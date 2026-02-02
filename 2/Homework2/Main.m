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

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[Femregion] = CreateFemregion(Data,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES
%==========================================================================

[D_no_bc,M_no_bc] = Matrix1D(Data,Femregion);

D = D_no_bc;
M = M_no_bc;
% BC
if strcmp(Data.boundary,'PP')

    % ---- Apply periodic constraint: dof(1) == dof(end)
    % Merge end DOF into first DOF (rows and columns)
    D(1,:) = D(1,:) + D(end,:);
    D(:,1) = D(:,1) + D(:,end);

    M(1,:) = M(1,:) + M(end,:);
    M(:,1) = M(:,1) + M(:,end);

    % Remove the redundant last DOF
    D = D(1:end-1, 1:end-1);
    M = M(1:end-1, 1:end-1);

elseif strcmp(Data.boundary,'RR')
    R = 0*D_no_bc;
    R(1,1) = Data.a;
    R(end,end) = Data.b;
    A_no_bc = D_no_bc + R;

elseif strcmp(Data.boundary,'AA')
    S = 0*D_no_bc;
    S(1,1) = Data.mu/Data.c * Data.alfa;
    S(end,end) = Data.mu/Data.c * Data.alfa;
end

%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================

[b_nbc] = Rhs1D(Data,Femregion);
b = b_nbc;

% --- Periodic reduction for RHS, consistent with M and D reduction
if strcmp(Data.boundary,'PP')
    b(1,:) = b(1,:) + b(end,:);
    b = b(1:end-1,:);
end

% ============================================================
% Build the full block matrices
% ============================================================
N  = size(M,1);
Z  = sparse(N,N);

Mfull = [M  Z;
         Z  M];

Sfull = [Z      -D;
        -Data.g*Data.H*D   Z];

% Forcing: usually only in first equation; set second block to zero
Bfull = [b; sparse(N, size(b,2))];

% ============================================================
% Add friction: q_t + gH eta_x = -gamma q
% => add +gamma*M*q in the second equation
% ============================================================
if isfield(Data,'gamma') && Data.gamma > 0
    Sfull = Sfull + [Z, Z;
                     Z, Data.gamma * M];
end

%==========================================================================
% BUILD INITIAL CONDITIONS AND TIME INTEGRATION (theta-method)
%==========================================================================

theta = 1/2;   % suggested

time = 0:Data.dt:Data.T;
nT   = numel(time);

% --- Build initial condition vectors eta0 and q0 at DOFs
x = Femregion.coord(:,1);   % adjust if your coord storage differs

eta0 = Data.eta0(x);
q0   = Data.q0(x);

% Periodic reduction for ICs (must match reduced M,D)
if strcmp(Data.boundary,'PP')
    eta0 = eta0(1:end-1);
    q0   = q0(1:end-1);
end

U = zeros(2*N, nT);
U(:,1) = [eta0; q0];

% --- Precompute matrices for theta-method:
% (M/dt + theta S) U_{k+1} = (M/dt - (1-theta) S) U_k - RHSpart
A = (1/Data.dt) * Mfull + theta * Sfull;
G = (1/Data.dt) * Mfull - (1-theta) * Sfull;

% --- Reflecting wall BC (q=0 at x=0,L): prepare modified system matrix
A_bc = A;
if strcmp(Data.boundary,'RW')
    idx_q_left  = N + 1;
    idx_q_right = 2*N;

    A_bc(idx_q_left, :) = 0;
    A_bc(idx_q_left, idx_q_left) = 1;

    A_bc(idx_q_right, :) = 0;
    A_bc(idx_q_right, idx_q_right) = 1;
end

% If you have forcing: use midpoint in time
% Here Bfull is (2N x nT). Use theta-weighted RHS:
% A U_{k+1} = G U_k + theta*B_{k+1} + (1-theta)*B_k
for k = 1:nT-1
    rhs = G * U(:,k) + theta * Bfull(:,k+1) + (1-theta) * Bfull(:,k);
    if strcmp(Data.boundary,'RW')
        rhs(idx_q_left)  = 0;
        rhs(idx_q_right) = 0;
        U(:,k+1) = A_bc \ rhs;
    else
        U(:,k+1) = A \ rhs;
    end
end

% Split solution
eta_snap = U(1:N, :);
q_snap   = U(N+1:end, :);

%==========================================================================
% POST-PROCESSING
%==========================================================================

% Movie / snapshots
if Data.snapshot
    t = 0;
    for i = 1:size(eta_snap,2)
        Snapshot(Femregion, eta_snap(:,i), q_snap(:,i), Data, t);
        % pause(0.01);
        t = t + Data.dt;
    end
end

% Solution at final time
Solutions = PostProcessing(Data, Femregion, eta_snap(:,end), q_snap(:,end));

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
    surf(tt, xx, eta_snap, 'EdgeColor', 'none');
    xlabel('time'); ylabel('x'); zlabel('\eta(x,t)');
    title('Surface plot: \eta');

    figure(4); clf;
    surf(tt, xx, q_snap, 'EdgeColor', 'none');
    xlabel('time'); ylabel('x'); zlabel('q(x,t)');
    title('Surface plot: q');
end

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
Errors = [];
if (Data.calc_errors)
    [Errors] = ComputeErrors(Data,Femregion,Solutions);
end



