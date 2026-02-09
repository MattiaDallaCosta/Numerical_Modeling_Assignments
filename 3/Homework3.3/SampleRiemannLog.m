function [h_s, u_s] = SampleRiemannLog(hs, us, hL, uL, hR, uR, g, h_ref)
% -------------------------------------------------------------------------
% Sample the exact Riemann solution at xi = x/t = 0 for the LOG-PRESSURE
% shallow water equations with P(h) = g*h*log(h/h_ref).
%
% Same logic as the standard sampler (check contact side, then wave type,
% then position relative to wave), but with:
%
%   - Celerity:     c(h) = sqrt(g*(1 + log(h/h_ref)))
%   - Shock speed:  S_K from general Hugoniot (not the sqrt(h*(h+hK)) form)
%   - Fan formula:  requires a mini-Newton solve (no closed form)
%
% INPUT:
%   - hs, us        Star-region state found by Newton
%   - hL, uL        Left cell state
%   - hR, uR        Right cell state
%   - g, h_ref      Gravity and reference depth
% OUTPUT:
%   - h_s, u_s      State at xi = 0
% -------------------------------------------------------------------------

    % Wave celerities
    wL = 1 + log(hL / h_ref);
    wR = 1 + log(hR / h_ref);
    ws = 1 + log(hs / h_ref);
    cL = sqrt(g * wL);
    cR = sqrt(g * wR);
    cs = sqrt(g * ws);

    if us >= 0
        % =================================================================
        % xi=0 is on the LEFT side of the contact --> check wave 1
        % =================================================================

        if hs <= hL
            % ----- LEFT RAREFACTION -----
            SHL = uL - cL;          % Head speed
            STL = us - cs;          % Tail speed

            if SHL >= 0
                % Entire fan to the right of xi=0
                h_s = hL;  u_s = uL;

            elseif STL <= 0
                % Entire fan to the left of xi=0 --> star region
                h_s = hs;  u_s = us;

            else
                % xi=0 is INSIDE the left rarefaction fan.
                %
                % At xi=0:  u - c = 0  -->  u = c = sqrt(g*w)
                % Riemann invariant (1-wave):
                %   u + (2*sqrt(g)/3)*w^(3/2) = uL + (2*sqrt(g)/3)*wL^(3/2)
                %
                % Substituting u = sqrt(g*w):
                %   sqrt(g*w) + (2*sqrt(g)/3)*w^(3/2) = RHS
                %
                % Solve for w by Newton iteration.

                RHS = uL + (2*sqrt(g)/3) * wL^(3/2);

                % Initial guess: average of wL and ws
                w = 0.5*(wL + ws);
                for iter = 1:20
                    sw = sqrt(w);
                    F  = sqrt(g)*sw*(1 + 2*w/3) - RHS;
                    dF = sqrt(g)*(1 + 2*w) / (2*sw);
                    dw = -F / dF;
                    w  = w + dw;
                    w  = max(w, 1e-12);
                    if abs(dw) < 1e-12, break; end
                end

                c_fan = sqrt(g * w);
                h_s   = h_ref * exp(w - 1);
                u_s   = c_fan;           % u = c at xi=0 for 1-wave
            end

        else
            % ----- LEFT SHOCK -----
            % General Hugoniot shock speed:
            %   SL = uL - sqrt[ (P(h*) - P(hL)) * h* / (hL*(h*-hL)) ]
            Ps = g * hs * log(hs / h_ref);
            PL = g * hL * log(hL / h_ref);
            SL = uL - sqrt((Ps - PL) * hs / (hL * (hs - hL)));

            if SL >= 0
                h_s = hL;  u_s = uL;
            else
                h_s = hs;  u_s = us;
            end
        end

    else
        % =================================================================
        % xi=0 is on the RIGHT side of the contact --> check wave 2
        % =================================================================

        if hs <= hR
            % ----- RIGHT RAREFACTION -----
            SHR = uR + cR;          % Head speed
            STR = us + cs;          % Tail speed

            if SHR <= 0
                % Entire fan to the left
                h_s = hR;  u_s = uR;

            elseif STR >= 0
                % Entire fan to the right --> star region
                h_s = hs;  u_s = us;

            else
                % xi=0 is INSIDE the right rarefaction fan.
                %
                % At xi=0:  u + c = 0  -->  u = -c = -sqrt(g*w)
                % Riemann invariant (2-wave):
                %   u - (2*sqrt(g)/3)*w^(3/2) = uR - (2*sqrt(g)/3)*wR^(3/2)
                %
                % Substituting u = -sqrt(g*w):
                %   sqrt(g*w) + (2*sqrt(g)/3)*w^(3/2) = (2*sqrt(g)/3)*wR^(3/2) - uR

                RHS = (2*sqrt(g)/3) * wR^(3/2) - uR;

                w = 0.5*(wR + ws);
                for iter = 1:20
                    sw = sqrt(w);
                    F  = sqrt(g)*sw*(1 + 2*w/3) - RHS;
                    dF = sqrt(g)*(1 + 2*w) / (2*sw);
                    dw = -F / dF;
                    w  = w + dw;
                    w  = max(w, 1e-12);
                    if abs(dw) < 1e-12, break; end
                end

                c_fan = sqrt(g * w);
                h_s   = h_ref * exp(w - 1);
                u_s   = -c_fan;          % u = -c at xi=0 for 2-wave
            end

        else
            % ----- RIGHT SHOCK -----
            Ps = g * hs * log(hs / h_ref);
            PR = g * hR * log(hR / h_ref);
            SR = uR + sqrt((Ps - PR) * hs / (hR * (hs - hR)));

            if SR <= 0
                h_s = hR;  u_s = uR;
            else
                h_s = hs;  u_s = us;
            end
        end
    end
end
