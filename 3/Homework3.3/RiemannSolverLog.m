function [fh, fq] = RiemannSolverLog(hL, uL, hR, uR, g, h_ref)
% -------------------------------------------------------------------------
% Exact Riemann solver for the LOG-PRESSURE shallow water equations:
%
%   F(U) = ( q ;  q^2/h + g*h*log(h/h_ref) )
%
% Same Newton iteration structure as the standard solver, but:
%   - Initial guess from two-rarefaction approximation with the new
%     Riemann invariant integral I(h) = (2*sqrt(g)/3)*(1+log(h/h_ref))^(3/2)
%   - State functions from StateFunctionsLog
%   - Sampling from SampleRiemannLog
%   - Final flux uses g*h*log(h/h_ref) instead of (1/2)*g*h^2
%
% INPUT:
%   - hL, uL    Left state
%   - hR, uR    Right state
%   - g          Gravity
%   - h_ref      Reference depth
% OUTPUT:
%   - fh, fq    Godunov flux components
% -------------------------------------------------------------------------

    tol = 1e-12;
    max_iter = 50;

    % Wave parameter w = 1 + log(h/h_ref)
    wL = 1 + log(hL / h_ref);
    wR = 1 + log(hR / h_ref);

    % Two-rarefaction initial guess:
    %   (4*sqrt(g)/3)*w*^(3/2) = (2*sqrt(g)/3)*(wL^(3/2)+wR^(3/2)) + (uL-uR)
    %   w* = [ (wL^(3/2)+wR^(3/2))/2 + 3*(uL-uR)/(4*sqrt(g)) ]^(2/3)
    %   h* = h_ref * exp(w* - 1)
    W = (wL^(3/2) + wR^(3/2))/2 + 3*(uL - uR)/(4*sqrt(g));
    W = max(W, 1e-6);     % ensure positivity
    ws = W^(2/3);
    hs = h_ref * exp(ws - 1);
    hs = max(hs, 1e-14);

    % Newton iteration on f_L(h*) + f_R(h*) + (uR - uL) = 0
    for iter = 1:max_iter
        [fL, dfL] = StateFunctionsLog(hs, hL, g, h_ref);
        [fR, dfR] = StateFunctionsLog(hs, hR, g, h_ref);

        res  = fL + fR + (uR - uL);
        dres = dfL + dfR;

        dh = -res / dres;
        hs = hs + dh;
        hs = max(hs, 1e-14);

        if abs(dh) / (abs(hs) + 1e-14) < tol, break; end
    end

    % Recover intermediate velocity
    [fL_val, ~] = StateFunctionsLog(hs, hL, g, h_ref);
    us = uL - fL_val;

    % Sample the Riemann solution at xi = 0
    [h_s, u_s] = SampleRiemannLog(hs, us, hL, uL, hR, uR, g, h_ref);

    % GODUNOV FLUX with log-pressure:
    %   F = ( h*u ;  h*u^2 + g*h*log(h/h_ref) )
    fh = h_s * u_s;
    fq = h_s * u_s^2 + g * h_s * log(h_s / h_ref);
end
