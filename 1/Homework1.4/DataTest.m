%% Data for the following PDE

function [Data] = DataTest(TestName)

if strcmp(TestName,'Homework1Point2')
    
    % Data name
    Data.name    = TestName;
    % Data domain
    Data.domain  = [0,1];

    % D = Dirichlet; N = Neumann; R = Robin; I = Impedance; P = periodic
    % We want: Neumann at x=0, Robin/impedance at x=1
    Data.boundary = 'ND';
    
    % Parameters and external forces 
    Data.ro    = 1;           % rho
    Data.mu    = 1;           % stiffness mu
    Data.omega = 5*pi;        % omega
    Data.vel   = 1;           % c = sqrt(mu/rho) = 1
    
    % Forcing term: f(x) = 0
    Data.force = @(x,omega,ro,vel) 0.*x;
    
    % Neumann BC at x = 0: u_x(0) = gN = omega
    Data.gN1   = @(omega,ro,vel) omega;
    
    % No Neumann BC at x = 1 (right BC is Robin/Impedance)
    Data.gN2   = @(omega,ro,vel) 0;
    
    % Robin / Impedance BC at x = 1:
    Data.gR1   = @(omega,ro,vel) 0;        % no Robin at x=0
    Data.gR2   = @(omega,ro,vel) 0;       % gA at x=1
    
    % (Dirichlet data not used in this test)
    Data.gD1   = @(omega,ro,vel) 0;
    Data.gD2   = @(omega,ro,vel) 0;
    
    % Exact solution for error analysis: u_ex(x) = sin(omega*x)
    Data.uex     = @(x,omega,ro,vel) sin(omega.*x);
    Data.graduex = @(x,omega,ro,vel) omega.*cos(omega.*x);
end

           



