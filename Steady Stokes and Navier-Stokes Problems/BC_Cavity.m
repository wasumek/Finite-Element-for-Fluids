function [bc, mesh] =  BC_Cavity(mesh, parameters)
% [Accd, bccd, velo] =  BC_Cavity(nx,ny,nunk)
% This function creates matrices Accd and bccd to impose Dirichlet 
% boudary conditions for the Cavity flow problem using Lagrange 
% multipliers method.
% Velocity is prescibed in the whole boundary and vertical component
% is zero everywhere whereas the horizontal one is zero everywhere 
% but the top, where it is equal one. 
% Besides, the function provides a velocity field velo which verifies boundary
% conditions and can be used as an initial guess for a nonlinear analysis.
%
% Input:
%   nx,ny:  number of elements on each direction
%   nunk:   number of degrees of freedom for the velocity field
%

% Imposed boundary conditions

D1u = ones(1,size(mesh.B1,2)) * parameters.u1(1);
D1v = ones(1,size(mesh.B1,2)) * parameters.u1(2);

D2u = ones(1,size(mesh.B2,2)) * parameters.u2(1);
D2v = ones(1,size(mesh.B2,2)) * parameters.u2(2);

D3u = ones(1,size(mesh.B3,2)) * parameters.u3(1);
D3v = ones(1,size(mesh.B3,2)) * parameters.u3(2);

D4u = ones(1,size(mesh.B4,2)) * parameters.u4(1);
D4v = ones(1,size(mesh.B4,2)) * parameters.u4(2);

Cy0 = [reshape([2*mesh.B1-1;    2*mesh.B1]  ,2*mesh.bc_n_nodes(1),1), ...
       reshape([D1u;            D1v      ]	,2*mesh.bc_n_nodes(1),1)];
Cx1 = [reshape([2*mesh.B2-1; 	2*mesh.B2] 	,2*mesh.bc_n_nodes(2),1), ...
       reshape([D2u;            D2v      ]  ,2*mesh.bc_n_nodes(2),1)];
Cx0 = [reshape([2*mesh.B4-1;  	2*mesh.B4] 	,2*mesh.bc_n_nodes(4),1), ...
       reshape([D4u;            D4v      ]  ,2*mesh.bc_n_nodes(4),1)];
Cy1 = [reshape([2*mesh.B3-1;   	2*mesh.B3] 	,2*mesh.bc_n_nodes(3),1), ...
       reshape([D3u;            D3v      ]  ,2*mesh.bc_n_nodes(3),1)];
C = [Cy0; Cx1; Cx0; Cy1];
bc.C = C;

% Boundary conditions' matrix
bc.n_Dir_nodes = size(C,1);
bc.Accd = zeros(bc.n_Dir_nodes,mesh.dfP);
bc.Accd(:,C(:,1)) = eye(bc.n_Dir_nodes);
% Boundary conditions' vector
bc.bccd = C(:,2);


