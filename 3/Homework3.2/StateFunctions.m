function [f, df] = StateFunctions(hs, hK, g)
% -------------------------------------------------------------------------
% State function f_K(h*) and its derivative df_K/dh* 
%
% Two branches:
%   h* <= h_K (rarefaction): f = 2*(sqrt(g*h*) - sqrt(g*h_K))
%
%   h* >  h_K (shock): f = (h* - h_K) * sqrt( (g/2)*(1/h* + 1/h_K) )
%
% INPUT:
%   - hs            gues of solution h*
%   - hK            h at right of left of the interface
%   - g             gravity
% OUTPUT:
%   - f             evaluation of the state function in hs
%   - df            evaluation of the derivated state function in hs
% -------------------------------------------------------------------------

    if hs <= hK
        % Rarefaction wave
        cs = sqrt(g * hs);
        cK = sqrt(g * hK);
        f  = 2 * (cs - cK);
        df = sqrt(g / hs);
    else
        % Shock wave
        A  = sqrt(0.5 * g * (1/hs + 1/hK));
        f  = (hs - hK) * A;
        dA = -g / (4 * A * hs^2);
        df = A + (hs - hK) * dA;
    end
end
