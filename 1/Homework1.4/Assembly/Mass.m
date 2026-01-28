function [M_loc]=Mass(dphiq,w_1D,nln,BJ)
%% [M_loc]=Mass(dphiq,w_1D,nln,BJ)
%==========================================================================
% Build the local mass matrix for the term (uv)
%==========================================================================
%    called in Matrix1D.m
%
%    INPUT:
%          dphiq       : (array real) evaluation of the basis function on
%                        quadrature nodes
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : Jacobian of the map 
%
%    OUTPUT:
%          M_loc       :  (array real) Local mass matrix

M_loc=zeros(nln,nln);

for i=1:nln
    for j=1:nln
        for k=1:length(w_1D)
            Binv = 1./BJ;      % inverse
            Jdet = BJ;         % determinant 
            M_loc(i,j) = M_loc(i,j) + (Jdet.*w_1D(k)) .* dphiq(1,k,i).* dphiq(1,k,j);
        end
    end
end

% --- Derivation of the local mass matrix for linear elements -------------
%
% Code formula:
%   M_loc(i,j) += (Jdet * w_1D(k)) * dphiq(1,k,i) * dphiq(1,k,j)
%
% where:
%   Jdet = BJ = h          (element length)
%   dphiq(1,k,1) = phi1(ξ_k)
%   dphiq(1,k,2) = phi2(ξ_k)
%
% Reference linear basis on [0,1]:
%   phi1(ξ) = 1 - ξ
%   phi2(ξ) = ξ
%
% --- M_loc(1,1) ----------------------------------------------------------
%   M11 = sum_k h * w_k * phi1(ξ_k)^2
%       ≈ h * ∫_0^1 (1 - ξ)^2 dξ
%       = h * ∫_0^1 (1 - 2ξ + ξ^2) dξ
%       = h * [ ξ - ξ^2 + ξ^3/3 ]_0^1
%       = h * (1 - 1 + 1/3)
%       = h/3
%
% --- M_loc(1,2) and M_loc(2,1) ------------------------------------------
%   M12 = M21 = sum_k h * w_k * phi1(ξ_k) * phi2(ξ_k)
%       ≈ h * ∫_0^1 (1 - ξ) ξ dξ
%       = h * ∫_0^1 (ξ - ξ^2) dξ
%       = h * [ ξ^2/2 - ξ^3/3 ]_0^1
%       = h * (1/2 - 1/3)
%       = h/6
%
% --- M_loc(2,2) ----------------------------------------------------------
%   M22 = sum_k h * w_k * phi2(ξ_k)^2
%       ≈ h * ∫_0^1 ξ^2 dξ
%       = h * [ ξ^3/3 ]_0^1
%       = h/3
%
% --- Final local mass matrix ---------------------------------------------
%   M_loc = (h/6) * [ 2  1 ;
%                     1  2 ];
%
% With Simpson's rule on [0,1]:
%   ξ = [0, 0.5, 1],  w = [1/6, 4/6, 1/6],
% the discrete sum in the code reproduces exactly these integrals.



                                              
