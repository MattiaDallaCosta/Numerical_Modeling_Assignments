function [A,M]=Matrix1D(Data,femregion)
%% [A,M] = Matrix1D(Data,Femregion)
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
%          A           : (sparse(ndof,ndof) real) stiffnes matrix
%          M           : (sparse(ndof,ndof) real) mass matrix


fprintf('Assembling the matrices M and A... \n');


% connectivity infos
ndof         = femregion.ndof;         % degrees of freedom (number of nodes)
nln          = femregion.nln;          % local degrees of freedom (2)
ne           = femregion.ne;           % number of elements (nEl)
connectivity = femregion.connectivity; % connectivity matrix


% shape functions in ref element: 1-xi and xi.
[basis] = ShapeBasis; 

% quadrature nodes and weights for integrals. (3) -> Cavalieri Simpson
[nodes_1D, w_1D] = Quadrature(3);

% evaluation of shape bases on quadrature nodes
[Phi,GradPhi] = EvalShapeBasis(basis,nodes_1D);

% Assembly begin ...
A = sparse(ndof,ndof);  % Global Stiffness matrix
M = sparse(ndof,ndof);  % Global mass matrix

% I assamble the matrix by calculating nEl local matrixes
for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
    
    % Remember we defined the basis in the reference element [0,1] (xi),
    % indeed to go from the normal integral in a generical interval
    % [x_i;x_{i+1}] we need a transformation map. The integral must indeed
    % include also the Jacobian factor.
    % Physical domain: (0,L)
    % Reference domain: (0,1)

    [BJ, xq] = GetJacobian(femregion.coord(iglo,:), nodes_1D);
    % BJ        = Jacobian of the elemental map (= h?)
    % pphys_1D  = vertex coordinates in the physical domain 
   
    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%
    % Evaluate at the nodes the stiffness mu
    muq = Data.mu(xq);

    % Local stiffness matrix (phi',phi')
    [A_loc] = Stiffness(GradPhi, w_1D, nln, BJ, muq);

    % Assembly phase for stiffness matrix
    A(iglo,iglo) = A(iglo,iglo) + A_loc; 
    
    %=============================================================%
    % MASS MATRIX
    %=============================================================%
    
    % Local mass matrix (phi,phi)
    [M_loc] = Mass(Phi, w_1D, nln, BJ);

    % Assembly phase for mass matrix
    M(iglo,iglo) = M(iglo,iglo) + M_loc;
end




