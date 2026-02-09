function sigma = MCslope(u_left, u_center, u_right, Delta_x)
% -------------------------------------------------------------------------
% Compute limited slope for MUSCL reconstruction. Using forward difference
% (u_{i+1} - u_i) / Delta_x
%
% INPUT:
%   - u_left        component h or q at previous cell
%   - u_center      component h or q at current cell    
%   - u_right       component h or q at next cell     
%   - Delta_x       spatial mesh size
% OUTPUT:
%   - sigma         limiter
% -------------------------------------------------------------------------

    fwd = (u_right  - u_center) / Delta_x;
    bwd = (u_center - u_left)   / Delta_x;

    % MC limiter
    if bwd == 0
        sigma = 0;
    else
        r = fwd / bwd;
        phi = max(0, min([2*r, 2, (1+r)/2]));
        sigma = phi * bwd;
    end
end