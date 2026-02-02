function [D,M]=Matrix1D(Data,femregion)
%% [D,M] = Matrix1D(Data,Femregion)
%==========================================================================
% Assembly of the stiffness matrices A and rhs f
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          femregion   : (struct)  see CreateFemregion.m
%
%    OUTPUT:
%          D           : (sparse(ndof,ndof) real) weak derivative matrix
%          M           : (sparse(ndof,ndof) real) mass matrix


fprintf('Assembling the matrices M and A... \n');


% connectivity infos
ndof         = femregion.ndof;         % degrees of freedom
nln          = femregion.nln;          % local degrees of freedom
ne           = femregion.ne;           % number of elements
connectivity = femregion.connectivity; % connectivity matrix


% shape functions
% [basis] = ShapeBasis;
 
% quadrature nodes and weights for integrals
% [nodes_1D, w_1D] = Quadrature(2);
 
% evaluation of shape bases on quadrature nodes
% [Phi,GradPhi] = EvalShapeBasis(basis,nodes_1D);

% Returns LGL nodes and weights on the reference element
[nodes_1D,w_1D] = xwlgl(nln);

% Evaluation of shape bases and their derivative on quadrature nodes
[Phi,GradPhi] = basis_and_der_at_lgl(nodes_1D,nln);



% Assembly begin ...
D = sparse(ndof,ndof);  % Global Stiffness matrix
M = sparse(ndof,ndof);  % Global mass matrix

for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
    
    % Compute the Jacobian
    [BJ, ~] = GetJacobian(femregion.coord(iglo,:), nodes_1D);
    % BJ        = Jacobian of the elemental map 
    % pphys_1D  = vertex coordinates in the physical domain 
   
    %=============================================================%
    % MASS MATRIX
    %=============================================================%
    % M_ij = ∫ phi_i * phi_j dx
    M_loc = Mass(Phi, w_1D, nln, BJ);
    M(iglo, iglo) = M(iglo, iglo) + M_loc;

    %=============================================================%
    % WEAK DERIVATIVE MATRIX
    %=============================================================%
    % D_ij = ∫ phi_j * (phi_i)_x dx
    D_loc = D_matrix(Phi, GradPhi, w_1D, nln, BJ);
    D(iglo, iglo) = D(iglo, iglo) + D_loc;
end


