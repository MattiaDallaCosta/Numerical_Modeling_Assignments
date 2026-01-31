%% Data for the following PDE

function [Data] = DataTest(TestName)

if strcmp(TestName,'Homework1Point2')
    
    % Data name
    Data.name    = TestName;
    % Data domain
    Data.domain  = [0,1.5];

    % D = Dirichlet; N = Neumann; R = Robin; I = Impedance; P = periodic
    % We want: Neumann at x=0, Robin/impedance at x=1
    Data.boundary = 'NR';
    
    % Parameters and external forces 
    % Density rho
    Data.ro = @(x) ...
        (x <= 1).*1 + ...
        (x > 1).*(1 ./ ((1 + 1i)));
    
    % Stiffness mu
    Data.mu = @(x) ...
        (x <= 0.5).*4 + ...
        (x > 0.5 & x <= 1).*1 + ...
        (x > 1).*((1 + 1i));        
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
    Data.gR2   = @(omega,ro,vel) 0;       % gA at x=1
    
    % (Dirichlet data not used in this test)
    Data.gD1   = @(omega,ro,vel) 0;
    Data.gD2   = @(omega,ro,vel) 0;
    
    % Exact solution for error analysis: u_ex(x) = sin(omega*x)
    omega = 5*pi;

    k1 = 5*pi/2;
    k2 = 5*pi;
    k3 = (5*pi/2)*(1 - 1i);
    
    A1 = (3 + 4i)/(50*pi);
    B1 = 1/(10*pi);
    
    a  = sqrt(2)/(25*pi);
    A2 = a*(1 - 2i);
    B2 = a*(-2 - 1i);
    
    C  = -0.000250864770314551 - 0.000752594310943653*1i;
    
    Data.uex = @(x,omega,ro,vel) ...
        (x <= 0.5).*( A1.*cos(k1*x) + B1.*sin(k1*x) ) + ...
        (x > 0.5 & x <= 1).*( A2.*cos(k2*x) + B2.*sin(k2*x) ) + ...
        (x > 1).*( C.*exp(-1i*k3.*(x - 1.5)) );

    Data.uexOld = @(x,omega,ro,vel) (x <= 0.5) .* ( ...
        (3+4i)/(50*pi) .* cos((5*pi/2).*x) + ...
        (1/(10*pi))      .* sin((5*pi/2).*x) ) ...
        + (x > 0.5) .* ( ...
        sqrt(2)*(1-2i)/(25*pi) .* cos(5*pi*x) - ...
        sqrt(2)*(2+1i)/(25*pi) .* sin(5*pi*x) );
    
    Data.graduex = @(x,omega,ro,vel) 0.*x;
end

           



