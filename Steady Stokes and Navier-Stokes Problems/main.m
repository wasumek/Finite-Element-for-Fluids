clear call; close all; clc; 
%% Definig the compuational domain.
% The computational domain is defined here. 

nsd = 2;                    % Number of spatial dimensions (2D only)
x0 = 0.0;                   % Left domain boundary
x1 = 1.0;                   % Right domain boundary
y0 = 0.0;                   % Bottom domain boundary
y1 = 1.0;                   % Top domain boundary
nx_elements = 20;           % Number of elements (ne) in x-direction
ny_elements = 20;           % Number of elements (ne) in y-direction
parameters.sine = 1;        % Use sine-spaced mesh.
parameters.elem_type = 2;   % Velocity Pressure combinations
                            %   1 -> Q1Q1 element
                            %   2 -> Q2Q1 element

parameters.matrix_fig       = true; 	% Plot matrix 
parameters.uvp_fig          = true; 	% Plot u,v and p fields

parameters.pshift = 0;   	% Pressure field offset
parameters.n_iter = 50;   	% Number of solver iterations

parameters.pspg = 0;        % PSPG stabilization on/off (1/0) 
                            % Not supported for full NS-equations

parameters.stokes = 0;   	% Stokes equations on/off (1/0)

%% Generate Mesh

    mesh = meshgen( x0, x1, y0, y1, nx_elements, ny_elements, nsd, parameters);
  
%% Defining the problem parameters
% Note that in this problem the fluid density is assumed to be unity.

    % Fluid viscosity  
    parameters.visc = 1.0/400;   
       
    % Initial solution guess (u,v)
    parameters.u_init = [0 0];
    
    % Body force term Fx, Fy
    parameters.s = [0 0];
                    
%% Setting boundary conditions  

    % Boundary condition flags for velocity field (u,v)
    % (1 if boundary is set 0 if not).
    parameters.w1 = [1 1];
    parameters.w2 = [1 1];
    parameters.w3 = [1 1];
    parameters.w4 = [1 1];

    % The Dirichlet value at the specific boundary can be set next. 
    parameters.u1 = [0.0 0.0];
    parameters.u2 = [0.0 0.0];
    parameters.u3 = [1.0 0.0];
    parameters.u4 = [0.0 0.0];


%% Solve solution using an iterative loop  
   
    [bc, sol] = solve_increment(mesh, parameters);
    
%% Plot solution
% Here the solution is being plotted. 

    plot_solution(mesh, bc, sol, parameters);


 





