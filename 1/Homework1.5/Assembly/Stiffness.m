function [K_loc] = Stiffness(Grad,w_1D,nln,BJ, mu)
%% [K_loc] = Stiffness(Grad,w_1D,nln,BJ)
%==========================================================================
% Build the local stiffness matrix for the term u'*v'
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          Grad        : (array real) evaluation of the derivative on
%                        quadrature nodes
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : (array real) Jacobian of the map 
%
%    OUTPUT:
%          K_loc       :  (array real) Local stiffness matrix


K_loc = zeros(nln,nln); % 2x2 matrix

%% General implementation -- to be used with general finite element spaces
for i=1:nln
    for j=1:nln
        for k=1:length(w_1D)
            Binv = 1./BJ;    % inverse 1/h
            Jdet = BJ;       % determinant h 
            K_loc(i,j) = K_loc(i,j) + mu(k) *(Jdet.*w_1D(k)) .* ( (Grad(k,:,i) * Binv) * (Grad(k,:,j) * Binv )');    
        end
    end
end

% --- Correct derivation of the stiffness matrix entries -------------------
%
% Starting from the implemented formula:
%   K_loc(i,j) += (Jdet * w_1D(k)) * ( (Grad(k,:,i) * Binv) * (Grad(k,:,j) * Binv) )
%
% where:
%   Jdet = h        (element length)
%   Binv = 1/h
%   Grad(k,:,1) = -1   for linear basis
%   Grad(k,:,2) =  1
%
% --- K_loc(1,1) ------------------------------------------------------------
%   K11 = sum_k [ h * w_k * (-1 * 1/h) * (-1 * 1/h) ]
%       = sum_k [ h * w_k * (1/h^2) ]
%       = (1/h) * sum_k w_k
% Using Simpson rule on [0,1]: sum_k w_k = 1
%       -> K11 = 1/h
%
% --- K_loc(1,2) ------------------------------------------------------------
%   K12 = sum_k [ h * w_k * (-1 * 1/h) * (1 * 1/h) ]
%       = sum_k [ h * w_k * (-1/h^2) ]
%       = -(1/h) * sum_k w_k
%       -> K12 = -1/h
%
% --- K_loc(2,1) ------------------------------------------------------------
% Symmetry of the bilinear form:
%       -> K21 = -1/h
%
% --- K_loc(2,2) ------------------------------------------------------------
%   K22 = sum_k [ h * w_k * (1 * 1/h) * (1 * 1/h) ]
%       = (1/h) * sum_k w_k
%       -> K22 = 1/h
%
% --- Final local stiffness matrix -----------------------------------------
%   K_loc = (1/h) * [ 1  -1 ;
%                    -1   1 ];
%
% This matches the standard 1D FEM stiffness matrix for linear elements.



                                              
                                              

