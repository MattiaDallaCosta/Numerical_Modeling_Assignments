function [Errors,Solutions,Femregion,Data] = MainHMZ(Data, nEl)
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

% Matrix of the system:
K = Data.omega^2 * M_nbc - A;

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

%==========================================================================
% PLOT
%==========================================================================

% Solve:
u_h = K\b;

% Plot
if (Data.visual_graph)
    [u_h] = PlotSol(Femregion, Data, u_h);
end


%==========================================================================
% SAVE SOLUTIONS
%==========================================================================

x = Femregion.coord(:,1);
u_ex = Data.uex(x,Data.omega,Data.ro,Data.vel);
Solutions = struct('u_ex',u_ex,'u_h',u_h,'x',x);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
Errors = [];
if (Data.calc_errors)
    [Errors] = ComputeErrors(Data,Femregion,Solutions);
end
