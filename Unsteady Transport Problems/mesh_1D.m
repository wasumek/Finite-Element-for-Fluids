function [ mesh ] = mesh_1D( x0, x1, n_elements, nsd, parameters)
% This function will generate a 1D mesh using the following strategy 
%
% Boundaries     --> x0                                  x1
% Mesh           --> [---0---0---0---0---0---0---0---0---]
% Node numbering --> 1   2   3   4   5   6   7   8   9  10
% The input is:
% x0            left boundary
% x1            right boundary
% n_elements    number of elements
% nsd           number of spatial dimensions
% fig           plotting flag

%% Mesh Characteristics
mesh.x0 = x0;                   % Left domain boundary
mesh.x1 = x1;                   % Right domain boundary
mesh.n_elements = n_elements;   % Number of elements
mesh.n_elm_nodes = 2;           % Number of nodes per element
mesh.nsd = nsd;                 % Number of spatial dimensions

%% Derived quantities:

mesh.L = mesh.x1 - mesh.x0;           % Domain length
mesh.h = mesh.L / mesh.n_elements;    % Element length.
mesh.hx = mesh.h;                     % Element lenght x-dir.
mesh.hy = 0.0;                        % Element length y-dir. (0 for 1D)
mesh.n_nodes = mesh.n_elements + 1;   % Number of nodes

%% Coordinate array:
% For the 1D case the dimension of the coordinate array is 
% (n_nodes,nsd). The row number of this array represents the node-number
% while the row entries represent the coordinates of this array. E.g. the
% entries of the 5th row contain the coordinates of the 5th node in your
% mesh. 

% Allocate coordinate array:
mesh.coord = zeros(mesh.n_nodes,1);

% Assemble coordinate array:
for i = 1:mesh.n_nodes;
    mesh.coord(i,1) = mesh.x0 + mesh.hx * (i-1);
end

%% Connectivity array:
% For the 1D case the dimension of the connectivity array is 
% (n_elements,n_elm_nodes). The row number of this array represents the
% element number while the row-entries contains the node-
% numbers of the element that is considered.

% Allocate connectivity array:
mesh.conn = zeros(mesh.n_elements,mesh.n_elm_nodes);

% Allocate connectivity array:
for i = 1:mesh.n_elements;
    mesh.conn(i,1) = i;
    mesh.conn(i,2) = i+1;
end

%% Boundary array:
% The boundary array contains information about which of the nodes are
% boundary nodes. In the 2D case there are two boundaries (i.e. boundary 1 
% and 2). The boundary array has size (n_nodes,1) and each row represents a
% node in the mesh. The array contains the boundary number when a node is
% part of a boundary. If a node is not part of a boundary the array entry
% is zero. 

% Allocate array:
mesh.boundary = zeros(mesh.n_nodes,1);

% Assemble array:
for i = 1:mesh.n_nodes;
    if mesh.coord(i,1) == x0;
        mesh.boundary(i,1) = 1;
    end
    if mesh.coord(i,1) == x1;
        mesh.boundary(i,1) = 2;        
    end
end

%% Print mesh
    if parameters.fig == true
        figure;

        %axis square                             % axes equally spaced
        title('Mesh with boundary nodes')       % Print title

        % plot mesh 
        % When using patch the array needs to have dimension (n_nodes,2). 
        % Hence we temporarily append a zero vector here.
        coordinates = [mesh.coord zeros(mesh.n_nodes,1)];

        xlabel('x');
        patch('Faces', mesh.conn, ...
              'Vertices', coordinates, ...
              'Marker','o', ...
              'FaceColor','w');

        hold on

        % plot boundary nodes in red
        for i = 1:mesh.n_nodes
            if mesh.boundary(i,1) == 1;
                x = mesh.coord(i,1);
                y = 0.0;
                plot(x,y,'Marker','o','Color','red',...
                     'MarkerFaceColor','r','MarkerSize',5)
            elseif mesh.boundary(i,1) == 2;
                x = mesh.coord(i,1);
                y = 0.0;
                plot(x,y,'Marker','o','Color','green',...
                     'MarkerFaceColor','green','MarkerSize',5)
            end

        end
    end
    
end