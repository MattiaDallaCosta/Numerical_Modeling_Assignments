function [D_loc] = D_matrix(Phi,Grad,w_1D,nln,BJ)
%% [D_loc] = D_matrix(Phi,Grad,w_1D,nln,BJ)
%==========================================================================
% Build the local D matrix for the term  ∫ phi_j * (phi_i)_x  dx
% i.e.  D_ij^loc = ∫ phi_j(x) * d/dx(phi_i(x)) dx
%==========================================================================
%
% INPUT:
%   Phi   : (array real) evaluation of basis functions on quadrature nodes
%           size typically (nq, nln)   (or compatible)
%   Grad  : (array real) evaluation of reference derivatives on quad nodes
%           size typically (nq, 1, nln) in 1D (as in your current code)
%   w_1D  : (array real) quadrature weights (nq)
%   nln   : (integer) number of local unknowns
%   BJ    : (array real) Jacobian of the map (scalar or (nq)-vector)
%
% OUTPUT:
%   D_loc : (array real) Local D matrix
%==========================================================================

nq = length(w_1D);
D_loc = zeros(nln,nln);
% dphiq(1,k,i).* dphiq(1,k,j)
for i = 1:nln
    for j = 1:nln
        for k = 1:nq

            Binv = 1 ./ BJ;   % inverse Jacobian (1D)
            Jdet = BJ;        % determinant Jacobian (1D)

            % physical derivative of phi_i: dphi_i/dx
            dphi_i_dx = (Grad(k,:,i) * Binv);   % scalar in 1D

            % value of phi_j at quad node
            phi_j = Phi(1,k,j);

            % D_ij = ∫ phi_j * dphi_i/dx dx
            D_loc(i,j) = D_loc(i,j) + (Jdet * w_1D(k)) * (phi_j * dphi_i_dx);
        end
    end
end
end
