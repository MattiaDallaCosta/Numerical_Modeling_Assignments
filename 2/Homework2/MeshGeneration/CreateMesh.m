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


x0 = Data.domain(1);
xL = Data.domain(2);

% Number of nodes per element
if Data.p == 1          % Linear case: 2 nodes per element
    npdx = 2;
elseif Data.p == 2
    npdx = 3;
elseif Data.p == 3
    npdx = 4;
elseif Data.p == 4      % Quadratic case: 5 nodes per element
    npdx = 5;
elseif Data.p == 5
    npdx = 6;
else
    disp('case not implemented')
end

% Mesh size
MeshSize = (xL-x0)./nEl;

% For each element find candidate nodes
i = 0; 
for ie = 1 : nEl
    % Right and left point of each element
    xb_ie = x0 + ie*MeshSize;
    xa_ie = xb_ie - MeshSize;
    
    % returns Legendre–Gauss–Lobatto nodes mapped to that element interval
    [xp,~] = xwlgl(npdx,xa_ie,xb_ie); 

    % Each element contains end and inner points, adjent elements have
    % indeed 1 same nodes. For this duplicates are preset
    p(i+1:npdx+i) = xp;
    i = i + npdx ;
end

% Delite duplicates
p = uniquetol(p,1.e-5);

% Number of global vertices
nVert = size(p,2); 

% Build connectivity
if (Data.p == 1)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]']';
elseif (Data.p == 2)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]']';
elseif (Data.p == 3)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]' [4:npdx-1:nVert-npdx+4]']';
elseif (Data.p == 4)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]' [4:npdx-1:nVert-npdx+4]' [5:npdx-1:nVert-npdx+5]']';
elseif (Data.p == 5)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]' [4:npdx-1:nVert-npdx+4]' [5:npdx-1:nVert-npdx+5]' [6:npdx-1:nVert-npdx+6]']';
else
    disp('case not implemented')
end

% Mesh data structure
Region = struct('dim',1,...
               'domain',Data.domain,...
               'h',MeshSize,...
               'nvert',nVert,...
               'ne',nEl,...
               'coord',p',...
               'boundary_points',[x0,xL],...
               'connectivity',t);
           
           
