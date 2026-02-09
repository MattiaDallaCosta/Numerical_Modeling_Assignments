function [f, df] = StateFunctionsLog(hs, hK, g, h_ref)
% -------------------------------------------------------------------------
% State function f_K(h*) and its derivative df_K/dh* for the LOG-PRESSURE
% shallow water equations with momentum flux:
%
%   F_2 = q^2/h + g*h*log(h/h_ref)
%
% The "pressure" is P(h) = g*h*log(h/h_ref), with P'(h) = g*(1+log(h/h_ref)).
% The wave celerity is c(h) = sqrt(P'(h)) = sqrt(g*(1+log(h/h_ref))).
%
% Two branches:
%
%   h* <= h_K (rarefaction):
%     The Riemann invariant integral gives:
%       f_K(h*) = (2*sqrt(g)/3) * [ w*^(3/2) - wK^(3/2) ]
%     where w(h) = 1 + log(h/h_ref).
%
%   h* > h_K (shock):
%     From the general Rankine-Hugoniot conditions:
%       f_K(h*) = sqrt[ (P(h*) - P(hK)) * (h* - hK) / (h* * hK) ]
%
% INPUT:
%   - hs        Guess of h*
%   - hK        Left or right state h_K
%   - g         Gravity
%   - h_ref     Reference depth (must be > 0)
% OUTPUT:
%   - f         f_K(h*)
%   - df        df_K/dh*
% -------------------------------------------------------------------------

    if hs <= hK
        % ==== RAREFACTION ====
        % w = 1 + log(h/h_ref);  c = sqrt(g*w)
        % Riemann invariant integral: I(h) = (2*sqrt(g)/3) * w^(3/2)
        % f_K = I(h*) - I(hK)
        ws = 1 + log(hs / h_ref);
        wK = 1 + log(hK / h_ref);

        f  = (2*sqrt(g)/3) * (ws^(3/2) - wK^(3/2));

        % df/dh* = (2*sqrt(g)/3) * (3/2)*ws^(1/2) * (1/hs)
        %        = sqrt(g*ws) / hs  =  c(h*) / h*
        df = sqrt(g * ws) / hs;

    else
        % ==== SHOCK (Rankine-Hugoniot) ====
        % P(h) = g*h*log(h/h_ref)
        % f^2 = (Ps - PK)*(hs - hK) / (hs*hK)
        Ps = g * hs * log(hs / h_ref);
        PK = g * hK * log(hK / h_ref);
        dP = Ps - PK;
        dh = hs - hK;

        Q  = dP * dh / (hs * hK);
        f  = sqrt(Q);

        % Derivative via dQ/dh*:
        %   Q = dP*dh / (hs*hK)
        %   dQ/dh* = [P'(h*)*dh + dP]/(hs*hK) - dP*dh/(hs^2*hK)
        Ps_prime = g * (1 + log(hs / h_ref));     % P'(h*)
        dQ = (Ps_prime * dh + dP) / (hs * hK) ...
           - dP * dh / (hs^2 * hK);
        df = dQ / (2 * f);
    end
end
