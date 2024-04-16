function [mesh] = mesh_2D(x0, x1, y0, y1, nx_elements, ny_elements, nsd, parameters)
% This function will generate a 2D mesh using the following strategy 
%
%      x0                    x1
%    y1 13----14-----15-----16 y1
%       |      |      |      |
%       |   7  |   8  |   9  |
%       |      |      |      |
%       9-----10-----11-----12
%       |      |      |      |
%       |   4  |   5  |   6  |
%       |      |      |      |
%       5------6------7------8
%       |      |      |      |
%       |   1  |   2  |   3  |
%       |      |      |      |
%    y0 1------2------3------4 y0
%      x0                    x1
% The input is:
% x0            left boundary
% x1            right boundary
% y0            bottom boundary
% y1            top boundary
% nx_elements   number of elements in x-direction
% ny_elements   number of elements in y-direction
% nsd           number of spatial dimensions
% fig           plotting flag


%% Mesh Characteristics
mesh.x0 = x0;                   % Left domain boundary                 
mesh.x1 = x1;                   % Right domain boundary
mesh.y0 = y0;                   % Bottom domain boundary
mesh.y1 = y1;                   % Top domain boundary
mesh.nx_elements = nx_elements; % Number of elements in x-direction
mesh.ny_elements = ny_elements; % Number of elements in y-direction
mesh.n_v_elm_nodes = 4;         % Number of velocity element nodes
mesh.n_p_elm_nodes = 4;         % Number of pressure element nodes
mesh.nsd = nsd;                 % Number of spatial dimensions
mesh.v_elem_type = 1;          	% Velocity element type 
mesh.p_elem_type = 1;          	% pressure element type 

%% Derived quantities: 
mesh.Lx = mesh.x1 - mesh.x0;                % Domain length in x-direction
mesh.Ly = mesh.y1 - mesh.y0;                % Domain length in y-direction
mesh.hx = mesh.Lx/mesh.nx_elements;         % Element length in x-direction
mesh.hy = mesh.Ly/mesh.ny_elements;         % Element length in y-direction
mesh.h = max(mesh.hx,mesh.hy);              % Maximum element length
mesh.nx_p_nodes = mesh.nx_elements + 1;  	% Num. of pressure nodes x-dir.
mesh.ny_p_nodes = mesh.ny_elements + 1; 	% Num. of pressure nodes y-dir.
mesh.nx_v_nodes = mesh.nx_elements + 1;  	% Num. of velocity nodes x-dir.
mesh.ny_v_nodes = mesh.ny_elements + 1; 	% Num. of velocity nodes y-dir.
mesh.n_elements = mesh.nx_elements ...      % Total number of elements 
                * mesh.ny_elements;  
mesh.n_v_nodes = mesh.nx_v_nodes ...        % Total num. of velocity nodes 
                * mesh.ny_v_nodes; 
mesh.n_p_nodes = mesh.nx_p_nodes ...        % Total num. of pressure nodes 
                * mesh.ny_p_nodes; 
            
mesh.edfU = mesh.nsd * mesh.n_v_elm_nodes; 	% Num. of velocity DOF per elm.
mesh.edfP = mesh.n_p_elm_nodes;             % Num. of pressure DOF per elm.
mesh.dfU = mesh.nsd * mesh.n_v_nodes;       % Num. of global velocity DOF
mesh.dfP = mesh.n_p_nodes;                  % Num. of global pressure DOF
            
%% Coordinate array.
% For the 2D case the dimension of the coordinate array is 
% (n_nodes,nsd). The row number of this array represents the node-number
% while the row entries represent the coordinates of this array. E.g. the
% entries of the 5th row contain the coordinates of the 5th node in your
% mesh. 

% Allocate coordinate array:
mesh.Vcoord = zeros(mesh.n_v_nodes,2);
xs = linspace(x0,x1,mesh.nx_v_nodes)';
ys = linspace(y0,y1,mesh.ny_v_nodes)';

if parameters.sine == 1
xs(2:mesh.nx_elements)  = (mesh.x1-mesh.x0)*(tanh(5/2)+tanh(5*(xs(2:mesh.nx_elements) ...
                            -(mesh.x0+mesh.x1)/2)))/(2*tanh(5/2))+mesh.x0;
ys(2:mesh.ny_elements)  = (mesh.y1-mesh.y0)*(tanh(5/2)+tanh(5*(ys(2:mesh.ny_elements) ...
                            -(mesh.y0+mesh.y1)/2)))/(2*tanh(5/2))+mesh.y0;
end                        
% Assemble coordinate array:
counter = 1;
for i = 1:mesh.ny_v_nodes
    for j = 1:mesh.nx_v_nodes
        mesh.Vcoord(counter,1) = xs(j); % x-coordinates
        mesh.Vcoord(counter,2) = ys(i); % y-coordinates
        counter = counter + 1;
    end
end

mesh.Pcoord = mesh.Vcoord;

%% Connectivity array:
% For the 2D case the dimension of the connectivity array is 
% (n_elements,n_elm_nodes). The row number of this array represents the
% element number while the row-entries contains the node-
% numbers of the element that is considered.

% Allocate connectivity array:
mesh.Vconn = zeros(mesh.n_elements,mesh.n_v_elm_nodes);

% Assemble connectivity array:
for i = 1:mesh.ny_elements
    for j = 1:mesh.nx_elements
        mesh.Vconn((i-1)*(nx_elements)+j,1)= (i-1)*mesh.nx_v_nodes + j;        % Bottom left node
        mesh.Vconn((i-1)*(nx_elements)+j,2)= (i-1)*mesh.nx_v_nodes + j + 1;    % Bottom right node
        mesh.Vconn((i-1)*(nx_elements)+j,3)= (i  )*mesh.nx_v_nodes + j + 1;    % Top right node
        mesh.Vconn((i-1)*(nx_elements)+j,4)= (i  )*mesh.nx_v_nodes + j;        % Top left node
    end
end

mesh.Pconn = mesh.Vconn;

%% Boundary data

% Allocate array:
mesh.bc_n_nodes = zeros(4,1);
counter = 1;

mesh.B1 = [];
mesh.B2 = [];
mesh.B3 = [];
mesh.B4 = [];

% Assemble array:
for i = 1:mesh.n_v_nodes
  	if mesh.Vcoord(i,2) == y0
        mesh.boundary(counter,1) = 1;
        mesh.boundary(counter,2) = i;
        mesh.B1 = [mesh.B1 i];
        mesh.bc_n_nodes(1) = mesh.bc_n_nodes(1) + 1;
        counter = counter + 1;
   	elseif mesh.Vcoord(i,2) == y1
        mesh.boundary(counter,1) = 3;
        mesh.boundary(counter,2) = i;
        mesh.B3 = [mesh.B3 i];
        mesh.bc_n_nodes(3) = mesh.bc_n_nodes(3) + 1;
        counter = counter + 1;
    elseif mesh.Vcoord(i,1) == x0
        mesh.boundary(counter,1) = 4;
        mesh.boundary(counter,2) = i;
        mesh.B4 = [mesh.B4 i];
        mesh.bc_n_nodes(4) = mesh.bc_n_nodes(4) + 1;
        counter = counter + 1;
	elseif mesh.Vcoord(i,1) == x1 
        mesh.boundary(counter,1) = 2;
        mesh.boundary(counter,2) = i;
        mesh.B2 = [mesh.B2 i];
        mesh.bc_n_nodes(2) = mesh.bc_n_nodes(2) + 1;
        counter = counter + 1;
    end

end

% Number of Dirichlet boundary nodes
mesh.n_Dir_nodes = counter - 1;


% Compute element length
hx  = xs(2:end) - xs(1:end-1);
hy  = ys(2:end) - ys(1:end-1); 
counter = 1;
for i = 1:mesh.nx_elements
    for j = 1:mesh.ny_elements
        mesh.h_elm(counter) = max(hx(i),hy(j));
        counter = counter + 1;
    end
end    

end