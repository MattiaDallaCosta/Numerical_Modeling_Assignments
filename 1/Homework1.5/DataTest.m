%% Data for the following PDE

function [Data] = DataTest(TestName)

if strcmp(TestName,'Homework1Point2')
    
    % Data name
    Data.name    = TestName;
    % Data domain
    Data.domain  = [0,1];

    % D = Dirichlet; N = Neumann; R = Robin; I = Impedance; P = periodic
    % We want: Neumann at x=0, Robin/impedance at x=1
    Data.boundary = 'NR';
    
    % Parameters and external forces 
    % Density rho
    Data.ro = @(x) ...
        (x <= 0.5).*1 + ...
        (x > 0.5 & x <= 1).*(1 ./ ((1 + 1i)));
    
    % Stiffness mu
    Data.mu = @(x) ...
        (x <= 0.5).*4 + ...
        (x > 0.5 & x <= 1).*(4*(1 + 1i));        
    Data.omega = 5*pi;        % omega
    Data.vel   = 1;           % c = sqrt(mu/rho) = 1
    
    % Forcing term: f(x) = 0
    Data.force = @(x,omega,ro,vel) 0.*x;
    
    % Neumann BC at x = 0: u_x(0) = gN = omega
    Data.gN1   = @(omega,ro,vel) 1;
    
    % No Neumann BC at x = 1 (right BC is Robin/Impedance)
    Data.gN2   = @(omega,ro,vel) 0;
    
    % Robin / Impedance BC at x = 1:
    Data.gR1   = @(omega,ro,vel) 0;        % no Robin at x=0
    Data.gR2   = @(omega,ro,vel) 0;        % gA at x=1
    
    % (Dirichlet data not used in this test)
    Data.gD1   = @(omega,ro,vel) 0;
    Data.gD2   = @(omega,ro,vel) 0;
    
    % Exact solution for error analysis:
    % Parameters
    % Parameters
    Data.omega = 5*pi;
    % Exact solution for the PML test (matches the PDE/BCs in the figure)

    x0 = 0.5;

    % Wavenumbers from (mu u')' + omega^2 rho u = 0  ->  k = omega*sqrt(rho/mu)
    k1 = Data.omega/2;                 % rho1/mu1 = 1/4
    k2 = Data.omega*(1-1i)/4;          % rho2/mu2 = 1/(4*(1+i)^2) = -i/8

    % Left region coefficients from mu(0) u_x(0) = 1 with mu1 = 4
    % u1(x) = A1 cos(k1 x) + B1 sin(k1 x),  u1'(0) = B1*k1
    B1 = 1/(4*k1);                % = 1/(2*omega)
    A1 = 1i*B1;                   % chosen to satisfy matching + outgoing condition

    % Enforce continuity at x0 and outgoing-wave condition in the PML layer
    u1_x0 = A1*cos(k1*x0) + B1*sin(k1*x0);

    % In PML: outgoing solution => u2(x) = A2*exp(-i*k2*x)
    % which is equivalent to A2 cos(k2 x) + B2 sin(k2 x) with B2 = -i A2
    A2 = u1_x0 * exp(1i*k2*x0);
    B2 = -1i*A2;

    % Piecewise exact solution
    Data.uex = @(x,omega,ro,vel) (x <= x0).*(A1.*cos(k1.*x) + B1.*sin(k1.*x)) + ...
        (x >  x0).*(A2.*cos(k2.*x) + B2.*sin(k2.*x));
end

           



