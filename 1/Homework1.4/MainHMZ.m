function [dispersion,h_lambda] = MainHMZ(Data, nEl, k_exact, nmodes)
%%
%    INPUT:
%          Data    : (struct) Data struct
%          nEl     : (int)    Number of mesh elements
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) finite element space
%
%          Data        : (struct)  Data struct
%

fprintf('============================================================\n')
fprintf(['Solving test ', Data.name, ' with ',num2str(nEl),' elements \n']);

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = CreateMesh(Data,nEl);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[Femregion] = CreateFemregion(Data,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES and RIGHT-HAND SIDE
%==========================================================================

[A_nbc,M_nbc] = Matrix1D(Data,Femregion);


%==========================================================================
% BUILD FINITE ELEMENTS RHS
%==========================================================================

[b_nbc] = Rhs1D(Data,Femregion);

%==========================================================================
% SOLVE LINEAR SYSTEM - COMPUTE BOUNDARY CONDITIONS
%==========================================================================

% Apply BC to b and A:
[A,b] = BoundaryConditions(A_nbc,b_nbc,Femregion,Data);
M_nbc(end,end) =1;
M_nbc(end,end-1) =0;
M_nbc(end-1,end) =0;
A(end,end) = 0;
A(end-1,end) = 0;
A(end,end-1) = 0;
% M_nbc = M_nbc(1:nEl,1:nEl);
% A = A(1:nEl,1:nEl);
if (Data.plot_matrixes)
    % Print the matrices
    fprintf('\n===============================\n');
    fprintf('  MASS MATRIX  M_{nbc}\n');
    fprintf('  (sparse display)\n');
    fprintf('===============================\n');
    disp(M_nbc);
    disp(full(M_nbc));
    
    fprintf('\n===============================\n');
    fprintf('  STIFFNESS MATRIX  A\n');
    fprintf('  (sparse display)\n');
    fprintf('===============================\n');
    disp(A);
    disp(full(A));
    
    fprintf('\n=====================\n');
    fprintf('      RHS vector b\n');
    fprintf('=====================\n');
    disp(b);
    disp(full(b));
end


%==========================================================================
% DISPERSION DATA
%==========================================================================
[V,D] = eigs(A, M_nbc, nmodes, 'sm');  % smallest eigenvalues

lambda_h = diag(D);
omega_h  = sqrt(lambda_h);
k_h = omega_h / Data.vel;

h_lambda = (Data.boundary(1)/nEl) * k_exact / (2*pi);
dispersion = k_h ./ k_exact;