%% Data for linearized shallow water equation
% EQ1: eta_{t}(x,t) - q_{x}(x,t) = 0  in (0,L) x (0,T] 
% EQ1: q_{t}(x,t) - gH * eta_{x}(x,t) = 0  in (0,L) x (0,T]
% IC:
%       eta(x,0) = e^(-50(x-0.5)^2);
%       q(x,0) = c * eta(x,0);
% BC: q(0,t) = q(L,t) and eta(0,t) = eta(L,t)

function [Data] = DataTest(TestName)

if strcmp(TestName,'Homework2.4')
    Data.name = TestName;

    % Domain (0,L)
    Data.L = 1;       
    Data.domain = [0, Data.L];

    % Boundary conditions: periodic in x
    % P = periodic;
    Data.boundary = 'PP';

    % Final time and time step
    Data.T  = 0.5;
    Data.dt = 0.0001;

    % Parameters
    Data.g  = 9.81;
    Data.H  = 1.0;

    % Wave speed: c = sqrt(gH)
    Data.c = sqrt(Data.g * Data.H);

    % No forcing
    Data.force = @(x,t) 0.*x.*t;

    % Initial conditions

    Data.eta0 = @(x) exp(-50*(x-0.5).^2);
    Data.gradeta0  = @(x) (-100*(x-0.5)).*exp(-50*(x-0.5).^2);

    Data.q0   = @(x) Data.c * Data.eta0(x);

    % Exact solution: eta(x,t) = eta(x - c*t, 0). I take the initial
    % condition and substitute x with x-c*t. mod(x - Data.c*t, Data.L)
    % permits to begin again from x = 0 when x-c*t > L.

    Data.etaex = @(x,t) Data.eta0( mod(x - Data.c*t, Data.L) );
    Data.gradetaex = @(x,t) Data.gradeta0(mod(x - Data.c*t, Data.L));

    Data.qex     = @(x,t) Data.c * Data.etaex(x,t);
    Data.gradqex = @(x,t) Data.c * Data.gradetaex(x,t);

elseif strcmp(TestName,'Homework2.5_1')
    %% Reflecting wall boundary conditions: q=0 at x=0,L

    Data.name = TestName;

    % Domain (0,L)
    Data.L = 1;
    Data.domain = [0, Data.L];

    % Boundary conditions for reflective wall
    Data.boundary = 'RW';

    % Final time and time step
    Data.T  = 0.5;
    Data.dt = 0.0001;

    % Polynomial degree
    Data.p = 4;

    % Parameters
    Data.g = 9.81;
    Data.H = 1.0;
    Data.c = sqrt(Data.g * Data.H);

    % No forcing
    Data.force = @(x,t) 0.*x;

    % Initial conditions (same pulse, right-going initially)
    Data.eta0 = @(x) exp(-50*(x-0.5).^2);
    Data.q0   = @(x) 0.*x;

    % No simple traveling-wave exact solution with reflecting walls:
    Data.etaex     = [];
    Data.gradetaex = [];
    Data.qex       = [];
    Data.gradqex   = [];

    % For BC convenience (if you keep these fields)
    Data.gD1 = @(t) 0.*t;  Data.gD2 = @(t) 0.*t;
    Data.gN1 = @(t) 0.*t;  Data.gN2 = @(t) 0.*t;
    Data.gR1 = @(t) 0.*t;  Data.gR2 = @(t) 0.*t;
    Data.gI1 = @(t) 0.*t;  Data.gI2 = @(t) 0.*t;

elseif strcmp(TestName,'Homework2.5_2')
    %% Input Data, periodic wave equation: eta_tt = gH eta_xx
    
    Data.name = TestName;

    % Friction
    Data.gamma = 1;

    % Domain (0,L)
    Data.L = 1;       
    Data.domain = [0, Data.L];

    % Boundary conditions: periodic in x
    % D = Dirichlet; P = periodic;
    Data.boundary = 'PP';
    %Data.boundary = 'RW';

    % Final time and time step
    Data.T  = 0.5;
    Data.dt = 0.0001;

    % --- Parameters to match eta_tt = gH eta_xx ---
    Data.g  = 9.81;
    Data.H  = 1.0;

    % wave speed: c = sqrt(gH)
    Data.c = sqrt(Data.g * Data.H);

    % No forcing
    Data.force = @(x,t) 0.*x.*t;

    % --- Choose a periodic exact solution for verification ---
    % Use a Fourier mode: eta(x,t) = cos(k x - w t),
    % periodic if k = 2*pi*m / L.
    % m = 1;                           % mode number (integer)
    % Data.k = 2*pi*m / Data.L;        % spatial wavenumber
    % Data.w = Data.c * Data.k;        % dispersion: w = c k

    % Initial conditions
    Data.eta0 = @(x) exp(-50*(x-0.5).^2);
    Data.q0   = @(x) Data.c * Data.eta0(x);
    %Data.q0   = @(x) 0.*x;

    % Exact solution requested: eta(x,t) = eta(x - c t, 0)
    % For periodic domain, wrap x-ct into [0,L)
    Data.etaex = @(x,t) Data.eta0( mod(x - Data.c*t, Data.L) );
    Data.gradeta0  = @(x) (-100*(x-0.5)).*exp(-50*(x-0.5).^2);
    Data.gradetaex = @(x,t) Data.gradeta0(mod(x - Data.c*t, Data.L));


    % (Optional) exact q if you want it consistent with q=c*eta traveling right
    Data.qex     = @(x,t) Data.c * Data.etaex(x,t);
    Data.gradqex = @(x,t) Data.c * Data.gradetaex(x,t);
end