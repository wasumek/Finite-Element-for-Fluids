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
mesh.n_elm_nodes = 4;           % Number of element nodes
mesh.nsd = nsd;                 % Number of spatial dimensions

%% Derived quantities: 
mesh.Lx = mesh.x1 - mesh.x0;           	% Domain length in x direction
mesh.Ly = mesh.y1 - mesh.y0;           	% Domain length in y direction
mesh.hx = mesh.Lx/mesh.nx_elements;   	% Element length in x direction
mesh.hy = mesh.Ly/mesh.ny_elements;   	% Element length in y direction
mesh.h = max(mesh.hx,mesh.hy);          % Maximum element length
mesh.nx_nodes = mesh.nx_elements + 1;  	% Number of nodes in x direction
mesh.ny_nodes = mesh.ny_elements + 1; 	% Number of nodes in y direction
mesh.n_elements = mesh.nx_elements ... 	% Total number of elements 
                * mesh.ny_elements;  
mesh.n_nodes = mesh.nx_nodes ...        % Total number of nodes 
                * mesh.ny_nodes;        
            
%% Coordinate array.
% For the 2D case the dimension of the coordinate array is 
% (n_nodes,nsd). The row number of this array represents the node-number
% while the row entries represent the coordinates of this array. E.g. the
% entries of the 5th row contain the coordinates of the 5th node in your
% mesh. 

% Allocate coordinate array:
mesh.coord = zeros(mesh.n_nodes,2);

% Assemble coordinate array:
counter = 1;
for i = 1:mesh.ny_nodes;
    for j = 1:mesh.nx_nodes;
        mesh.coord(counter,1) = mesh.y0 + (j-1) * mesh.hx; % x-coordinates
        mesh.coord(counter,2) = mesh.x0 + (i-1) * mesh.hy; % y-coordinates
        counter = counter + 1;
    end
end

%% Connectivity array:
% For the 2D case the dimension of the connectivity array is 
% (n_elements,n_elm_nodes). The row number of this array represents the
% element number while the row-entries contains the node-
% numbers of the element that is considered.

% Allocate connectivity array:
mesh.conn = zeros(mesh.n_elements,mesh.n_elm_nodes);

% Assemble connectivity array:
for i = 1:mesh.ny_elements;
    for j = 1:mesh.nx_elements;
        mesh.conn((i-1)*(nx_elements)+j,1)= (i-1)*mesh.nx_nodes + j;        % Bottom left node
        mesh.conn((i-1)*(nx_elements)+j,2)= (i-1)*mesh.nx_nodes + j + 1;    % Bottom right node
        mesh.conn((i-1)*(nx_elements)+j,3)= (i  )*mesh.nx_nodes + j + 1;    % Top right node
        mesh.conn((i-1)*(nx_elements)+j,4)= (i  )*mesh.nx_nodes + j;        % Top left node
    end
end

%% Boundary array
% The boundary array contains information about which of the nodes are
% boundary nodes. In the 2D case there are four boundaries (the bot, 
% right, top and left wall numbers 1 to 4 respectively). The boundary array
% has size (n_nodes,1) and each row represents a node in the mesh. The 
% array contains the boundary number when a node is part of a boundary. If 
% a node is not part of a boundary the array entry is zero. 

% Allocate array:
mesh.boundary = zeros(mesh.n_nodes,1);

% Assemble array:
for i = 1:mesh.n_nodes

    if mesh.coord(i,1) == x0;
        mesh.boundary(i,1) = 4;
    elseif mesh.coord(i,2) == y1;
        mesh.boundary(i,1) = 3;    
    elseif mesh.coord(i,1) == x1 
        mesh.boundary(i,1) = 2;
    elseif mesh.coord(i,2) == y0;
        mesh.boundary(i,1) = 1;
    end

end

%% Print mesh
if parameters.fig == true;
    figure

    axis square                             % axes equally spaced
    title('Mesh with boundary nodes')       % Print title

    % plot mesh 
    patch('Faces', mesh.conn, ...
          'Vertices', mesh.coord, ...
          'Marker','o', ...
          'FaceColor','w');

    hold on

    % plot boundary
    for i = 1:mesh.n_nodes;
        if mesh.boundary(i) == 1;
            x = mesh.coord(i,1);
            y = mesh.coord(i,2);
            plot(x,y,'Marker','o','Color','red',...
                 'MarkerFaceColor','r','MarkerSize',5);
        elseif mesh.boundary(i) == 2;
            x = mesh.coord(i,1);
            y = mesh.coord(i,2);
            plot(x,y,'Marker','o','Color','blue',...
                 'MarkerFaceColor','blue','MarkerSize',5);
        elseif mesh.boundary(i) == 3;
            x = mesh.coord(i,1);
            y = mesh.coord(i,2);
            plot(x,y,'Marker','o','Color','green',...
                 'MarkerFaceColor','green','MarkerSize',5);
        elseif mesh.boundary(i) == 4;
            x = mesh.coord(i,1);
            y = mesh.coord(i,2);
            plot(x,y,'Marker','o','Color','yellow',...
                 'MarkerFaceColor','yellow','MarkerSize',5);
        end
    end
end

end