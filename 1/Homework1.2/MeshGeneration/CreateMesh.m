function [Region] = CreateMesh(Data,nEl)
%% [Region] = C_create_mesh(Data, nEl)
%=========================================================================
% Creates regular mesh
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data    : (struct)  see DataTest.m
%          nEl     : (int)    Number of mesh elements  
%
%    OUTPUT:
%          Region  : (struct) having fields: dimension
%                                            domain 
%                                            mesh size
%                                            number of vertices
%                                            number of elements
%                                            coordinates
%                                            boundary points
%                                            connectivity


% Boundary edges
x0 = Data.domain(1);
xL = Data.domain(2);


%================================================
% Geometrical info
 nVert = nEl + 1;                   % # of points = # of interval +1
 p = linspace(x0,xL,nVert);         % equispaced points in the domain
 t = [[1:nVert-1]' [2:nVert]']';    % first column: [1 2 ...]
                                    % second column: [2 3 ...]
                                    % Connectivity: (12) (23) (34)
 MeshSize = (xL-x0)./nEl;           % equispaced points in the domain
%================================================

% Mesh data structure
Region = struct('dim',1,...
               'domain',Data.domain,...
               'h',MeshSize,...
               'nvert',nVert,...
               'ne',nEl,...
               'coord',p',...
               'boundary_points',[x0,xL],...
               'connectivity',t);
           
           
