function [h_s, u_s] = SampleRiemann(hs, us, hL, uL, hR, uR, g)
% -------------------------------------------------------------------------
% Sample the exact Riemann solution at xi = x/t = 0 to determine
% which wave region the interface x_{i+1/2} falls into.
%
% The Riemann solution has the following structure:
%
%   U_L ---[wave 1]--- U* ---[wave 2]--- U_R
%
% Newton's method has already found (hs, us) = (h*, u*), the state
% in the star region between the two waves. But the Godunov flux
% needs the state at xi = 0 (the cell interface), which may or may
% not be in the star region depending on the wave speeds.
%
% This function checks where xi = 0 falls and returns the correct state.
%
% INPUT:
%   - hs, us        Star-region state (h*, u*)
%   - hL, uL        Left cell state  (cell i)
%   - hR, uR        Right cell state (cell i+1)
%   - g             Gravitational acceleration
%
% OUTPUT:
%   - h_s, u_s      The actual state at xi = 0, which can be:
%                      (hL, uL)  if both waves went right
%                      (hs, us)  if the interface is in the star region
%                      (hR, uR)  if both waves went left
%                      (h_fan, u_fan) if xi=0 is inside a rarefaction fan
% -------------------------------------------------------------------------

    % First, compute the wave celerities for left, right, and star states.
    % These are needed to determine wave speeds.
    cL = sqrt(g * hL);      % Wave celerity of the left state
    cR = sqrt(g * hR);      % Wave celerity of the right state
    cs = sqrt(g * hs);      % Wave celerity of the star state

    % =====================================================================
    % STEP 1: Check the sign of u* to determine which SIDE of the contact
    %         discontinuity the interface xi=0 lies on.
    %
    %   u* >= 0  means the contact moves to the right,
    %            so xi=0 is on the LEFT side --> we only need to check wave 1
    %
    %   u* < 0   means the contact moves to the left,
    %            so xi=0 is on the RIGHT side --> we only need to check wave 2
    % =====================================================================

    if us >= 0
        % =================================================================
        % xi=0 is on the LEFT side of the contact.
        % We must now determine what wave 1 (the left wave) is doing.
        %
        % STEP 2a: Is wave 1 a rarefaction or a shock?
        %   h* <= hL  -->  rarefaction (water depth decreases across wave)
        %   h* >  hL  -->  shock       (water depth jumps up)
        % =================================================================

        if hs <= hL
            % =============================================================
            % CASE: LEFT RAREFACTION
            %
            % A rarefaction is not a single wavefront but a continuous fan
            % spread between a head (fastest edge) and a tail (slowest edge).
            %
            %   SHL = head speed = uL - cL  (leftmost edge of the fan)
            %   STL = tail speed = u* - c*  (rightmost edge of the fan)
            %
            % Three sub-cases depending on where xi=0 sits relative to the fan:
            %
            %        SHL        STL
            %         \   fan   \
            %          \ ~~~~~~ \
            %           \        \
            %   U_L      \  ???  \   U*
            %             \       \
            %    ----------●-----------  xi
            %              0
            % =============================================================
            SHL = uL - cL;       % Head speed (leftmost edge of fan)
            STL = us - cs;       % Tail speed (rightmost edge of fan)

            if SHL >= 0
                % ---------------------------------------------------------
                % Sub-case 1: The head of the fan is already to the right
                % of xi=0. The entire fan is to our right.
                % We haven't been reached by any wave yet.
                %
                %          SHL>0   STL
                %   U_L     \  fan  \   U*
                %    (us)    \       \
                %   ---●------\-------\---
                %      0
                %
                % Result: xi=0 sees the original left state.
                % ---------------------------------------------------------
                h_s = hL;
                u_s = uL;

            elseif STL <= 0
                % ---------------------------------------------------------
                % Sub-case 2: The tail of the fan has already passed to
                % the left of xi=0. The entire fan is to our left.
                % We are fully inside the star region.
                %
                %     SHL    STL<0
                %      \  fan  \    U*
                %       \       \  (us)
                %   -----\-------\---●---
                %                    0
                %
                % Result: xi=0 sees the star state.
                % ---------------------------------------------------------
                h_s = hs;
                u_s = us;

            else
                % ---------------------------------------------------------
                % Sub-case 3: SHL < 0 < STL
                % The head has passed us but the tail hasn't yet.
                % We are INSIDE the rarefaction fan.
                %
                %     SHL    0    STL
                %      \    ●     \
                %       \ (us) ~~~~\
                %        \~~~~~~~~~~\
                %
                % Inside the fan, the solution varies continuously.
                % For the left wave family (1-rarefaction), the Riemann
                % invariant u + 2c = const gives:
                %   u + 2c = uL + 2cL   (constant across the fan)
                % and the self-similar condition xi = u - c = 0 gives:
                %   u = c   (at xi=0)
                % Solving:  3c = uL + 2cL  -->  c_fan = (uL + 2cL)/3
                % ---------------------------------------------------------
                c_fan = (uL + 2*cL) / 3;
                h_s = c_fan^2 / g;
                u_s = c_fan;
            end

        else
            % =============================================================
            % CASE: LEFT SHOCK
            %
            % A shock is a single discontinuity travelling at speed SL,
            % computed from the Rankine-Hugoniot conditions.
            % Only two sub-cases: xi=0 is either left or right of it.
            %
            %          SL
            %          |
            %   U_L    |   U*
            %          |
            %   -------●-----------  xi
            %          0
            % =============================================================
            SL = uL - cL * sqrt(hs*(hs + hL) / (2*hL^2));

            if SL >= 0
                % ---------------------------------------------------------
                % The shock hasn't reached xi=0 yet (it's to our right).
                % We still see the original left state.
                % ---------------------------------------------------------
                h_s = hL;
                u_s = uL;
            else
                % ---------------------------------------------------------
                % The shock has already passed xi=0 (it's to our left).
                % We are in the star region.
                % ---------------------------------------------------------
                h_s = hs;
                u_s = us;
            end
        end

    else
        % =================================================================
        % u* < 0: xi=0 is on the RIGHT side of the contact.
        % Everything is symmetric: we check wave 2 (the right wave).
        %
        % STEP 2b: Is wave 2 a rarefaction or a shock?
        %   h* <= hR  -->  rarefaction
        %   h* >  hR  -->  shock
        % =================================================================

        if hs <= hR
            % =============================================================
            % CASE: RIGHT RAREFACTION
            %
            % Same logic as left rarefaction but mirrored.
            % Head speed is now the rightmost edge, tail is the leftmost.
            %
            %   SHR = head speed = uR + cR  (rightmost edge of fan)
            %   STR = tail speed = u* + c*  (leftmost edge of fan)
            % =============================================================
            SHR = uR + cR;       % Head speed (rightmost edge of fan)
            STR = us + cs;       % Tail speed (leftmost edge of fan)

            if SHR <= 0
                % ---------------------------------------------------------
                % Entire fan is to our left. We see the original right state.
                % ---------------------------------------------------------
                h_s = hR;
                u_s = uR;

            elseif STR >= 0
                % ---------------------------------------------------------
                % Entire fan is to our right. We are in the star region.
                % ---------------------------------------------------------
                h_s = hs;
                u_s = us;

            else
                % ---------------------------------------------------------
                % xi=0 is inside the right rarefaction fan.
                % For the right wave family (2-rarefaction):
                %   u - 2c = uR - 2cR   (Riemann invariant)
                %   u + c = 0            (self-similar at xi=0)
                % Solving:  3c = 2cR - uR  -->  c_fan = (2cR - uR)/3
                % And u = -c_fan (negative because right-going fan sampled
                % at xi=0 when the fan straddles the origin)
                % ---------------------------------------------------------
                c_fan = (2*cR - uR) / 3;
                h_s = c_fan^2 / g;
                u_s = -c_fan;
            end

        else
            % =============================================================
            % CASE: RIGHT SHOCK
            %
            % Single discontinuity at speed SR (Rankine-Hugoniot).
            % Symmetric to the left shock case.
            % =============================================================
            SR = uR + cR * sqrt(hs*(hs + hR) / (2*hR^2));

            if SR <= 0
                % ---------------------------------------------------------
                % Shock is to our left. We see the original right state.
                % ---------------------------------------------------------
                h_s = hR;
                u_s = uR;
            else
                % ---------------------------------------------------------
                % Shock is to our right. We are in the star region.
                % ---------------------------------------------------------
                h_s = hs;
                u_s = us;
            end
        end
    end
end