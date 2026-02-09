function [fh, fq] = RiemannSolver(hL, uL, hR, uR, g)
% -------------------------------------------------------------------------
% Exact Riemann solver for the 1D shallow water equations. Using Newton's 
% method.
%
% INPUT:
%   - hL            h at left of the interfacace
%   - qL            q at right of the interfacace
%   - hR            h at left of the interfacace
%   - qR            q at right of the interfacace
%   - g             gravity
% OUTPUT:
%   - fh            first component of Godunov flux
%   - fq            second component of Godunov flux 
% -------------------------------------------------------------------------


    tol = 1e-12;    % error tollerance
    max_iter = 50;  % max iteration

    % c_L and c_R
    cL = sqrt(g * hL);  
    cR = sqrt(g * hR);

    % Initial guess:
    hs = ((cL + cR + 0.5*(uL - uR)) / (2*sqrt(g)))^2;
    hs = max(hs, 1e-14); % positivity of h

    % Newton iteration on f_L(h*) + f_R(h*) + (u_R - u_L) = 0
    for iter = 1:max_iter
        % Evaluate state functions and their derivative
        [fL, dfL] = StateFunctions(hs, hL, g);
        [fR, dfR] = StateFunctions(hs, hR, g);
        
        % Error
        res  = fL + fR + (uR - uL);

        % Derivatives of the state functions
        dres = dfL + dfR;

        % Update the initial guess hs = hs - res/dres
        dh = -res / dres;
        hs = hs + dh;
        hs = max(hs, 1e-14);

        % Check if error is below treshold
        if abs(dh) / (abs(hs) + 1e-14) < tol
            break;
        end
    end

    % Recover intermediate velocity
    [fL_val, ~] = StateFunctions(hs, hL, g);
    us = uL - fL_val;

    % Sample the Riemann solution at x/t = 0 to get the Godunov state
    [h_s, u_s] = SampleRiemann(hs, us, hL, uL, hR, uR, g);

    % GODUNOV FLUX = physical flux evaluated at U*(0)
    fh = h_s * u_s;
    fq = h_s * u_s^2 + 0.5 * g * h_s^2;
end