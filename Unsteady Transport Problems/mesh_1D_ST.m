function [mesh] = mesh_1D_ST(x0, x1, nx_elements, parameters)
% This function will generate a 2D mesh using the following strategy 
%
%      x0                    x1
%       5------6------7------8 upper time level
%       |      |      |      |
%       |   1  |   2  |   3  |  single time slab
%       |      |      |      |
%       1------2------3------4 lower time level
%      x0                    x1
% The input is:
% x0            left boundary
% x1            right boundary
% y0            bottom boundary
% y1            top boundary
% nx_elements   number of elements in x-direction
% ny_elements   number of elements in y-direction
% nsd           number of spatial dimensions
% parameters    various numerical parameters


%% Mesh Characteristics
mesh.x0 = x0;                   % Left domain boundary                 
mesh.x1 = x1;                   % Right domain boundary
mesh.nx_elements = nx_elements; % Number of elements in space-time
mesh.n_elm_nodes_s = 2;     % Number of element nodes in space only
mesh.n_elm_nodes = 4;           % Number of element nodes
mesh.nsd = 1;                 % Number of spatial dimensions

%% Derived quantities: 
mesh.Lx = mesh.x1 - mesh.x0;           	% Domain length in x direction
mesh.hx = mesh.Lx/mesh.nx_elements;   	% Element length in x direction
mesh.hy = 0;                            % Element length in x direction
mesh.h = mesh.hx;                       % Maximum element length
mesh.nx_nodes = mesh.nx_elements + 1;  	% Number of nodes in x direction
mesh.n_elements = mesh.nx_elements; 	% Total number of nodes 

mesh.n_nodes = mesh.nx_nodes * 2;       % Total number of nodes 
                      
            
%% Coordinate array.

% Allocate array:
mesh.coord = zeros(mesh.n_nodes,mesh.nsd);

% Assemble array:
counter = 1;
for j = 1:mesh.nx_nodes;
    
    % x-coordinates lower time-slab
    mesh.coord(counter,1) = mesh.x0 + (j-1) * mesh.hx; 
   
    % x-doordinates at upper time-slab (simply an offset of the lower
    % time-slab coordinates).
    mesh.coord(counter + mesh.nx_nodes,1) = mesh.coord(counter,1);
    
    counter = counter + 1;
end

%% Connectivity array:

% Allocate array:
mesh.conn = zeros(mesh.n_elements,mesh.n_elm_nodes);

% Assemble array:
for i = 1:mesh.nx_elements;
    for j = 1:mesh.nx_elements;
        mesh.conn(j,1)= j;        % Bottom left node
        mesh.conn(j,2)= j+ 1;    % Bottom right node
        mesh.conn(j,3)= mesh.nx_nodes + j + 1;    % Top right node
        mesh.conn(j,4)= mesh.nx_nodes + j;        % Top left node
    end
end

%% Boundary array

% Allocate array:
mesh.boundary = zeros(mesh.n_nodes,1);

% Assemble array:
for i = 1:mesh.n_nodes

    if mesh.coord(i,1) == x0;
        mesh.boundary(i,1) = 1; 
    elseif mesh.coord(i,1) == x1 
        mesh.boundary(i,1) = 2;
    end

end

%% Time slab marker array

% This array contains a 1 for the nodes that are part of the lower
% time-slab and a 2 for nodes that are part of the upper time-slab.

% Allocate array:
mesh.time_slab = zeros(mesh.n_nodes,1);

% Assemble array:
for i = 1:mesh.n_nodes
    
    if i <= mesh.nx_nodes;
        mesh.time_slab(i,1) = 1; 
    else
        mesh.time_slab(i,1) = 2;
    end

end

%% Print mesh
    if parameters.fig == true
        figure;

        %axis square                             % axes equally spaced
        title('Space-only mesh with boundary nodes')       % Print title

        % plot mesh 
        % When using patch the array needs to have dimension (n_nodes,2). Hence
        % we temporarily append a zero vector here.
        coordinates = [mesh.coord(1:mesh.nx_nodes,1) zeros(mesh.nx_nodes,1)];

        xlabel('x');
        patch('Faces', mesh.conn(:,1:2), ...
              'Vertices', coordinates, ...
              'Marker','o', ...
              'FaceColor','w');

        hold on

        % plot boundary nodes in red
        for i = 1:mesh.nx_nodes
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