clear all; close all; clc; 
%% Definig the compuational domain.
% The computational domain is defined here. For the 1D case all the
% variables related to the y-direction are ignored.

%% Mesh data
nsd = 1;                    % Number of spatial dimensions
x0 = 0.0;                   % Left domain boundary
x1 = 1.0;                   % Right domain boundary
y0 = 0.0;                   % Bottom domain boundary
y1 = 1.0;                   % Top domain boundary
nx_elements = 150;           % Number of elements (ne) in x-direction
ny_elements = 10;           % Number of elements (ne) in x-direction
parameters.gauss_int = 2;   % For project 2 we will use 2 Gauss points per 
                            % element in 1D. This is set simply by: 
                            % parameters.gauss_int = 2;

%% Postprocessing
parameters.fig = false;     % Plot grid "true"/"false"
parameters.exact = 'fig516'; % For Problem 5.5 set this variable to
                            % parameters.exact = 'fig55'
                            % For Problem 5.14 set this variable to
                            % parameters.exact = 'fig514'
                            % For Guassian hill (5.16) set this variable to
                            % parameters.exact = 'fig516';

%% Time integration
parameters.space_time = 0;  % 1 --> Use space-time formulation (1D only) 
                            % 0 --> Use theta-method (Crank-Nicolson)
parameters.theta = 0.5;     % Theta-method selector.                            
parameters.dt = 0.02;     	% Time-step size.
parameters.tot_t_iter = 50*0.6; % Number of time steps



%% Generate Mesh

    mesh = meshgen( x0, x1, y0, y1, nx_elements, ny_elements, nsd, parameters);
        
    % Note that the mesh parameters initially given as imput arguments for
    % this routine are now stored in the "mesh" structure.

%% Select stabilization method

    % case 0 => No stabilization
    % case 1 => Full upwinding
    % case 2 => Balancing diffusion
    % case 3 => SUPG 
    % case 4 => GLS

    parameters.stab = 0;
        
%% Defining the problem parameters
% Here the diffusion coefficient and source term are specified and stored 
% in the structure "parameters"

    % Diffusion coefficient
    parameters.diffusion = 1/200/150; 
    %1/10=>Pe0.5, 1/100=>Pe5 for 10 elements
    %1/150=>Pe0.5, 1/1500=>Pe5 for 150 elements
    
    % Advection "flow" field (uniform over the domain)
    % For 1D case the second vector component must be set to zero.
    parameters.advection = [1; 0];  

    % Define source term
    parameters.s = 0.0;

    % Compute Peclet number Pe_x Pe_y and Pe_norm
    parameters = compute_peclet(mesh, parameters);

    % Compute Courant number
    parameters = compute_courant(mesh, parameters);
                    
%% Setting boundary and initial conditions  

    % Define initial condition
    parameters.u_init = 'fig516';    % Set this parameters to 
                                % parameters.u_init = 'fig516'; for initial
                                % condition of Gaussian hill problem. 
    
    % Boundary condition flags (1 if boundary is set 0 if not).
    parameters.w1 = 1;
    parameters.w2 = 1;
    parameters.w3 = 1;
    parameters.w4 = 1;

    % The value at of the specific boundary can be set next. 
    parameters.u1 = 0;
    parameters.u2 = 0;
    parameters.u3 = 0;
    parameters.u4 = 1;

%% Time solve

    % Timestep loop
    % For each timestep an FE problem is solved that uses either the
    % initial or current contition.
    sol = solve_time(mesh,parameters);

%% Print simulation data

	print_status(parameters, sol, 6);
 
%% Post processing solution

    % Maximum fps during plotting of solution
    parameters.fps = 10;
    
    % Plot the solution
    figure2 = plot_solution(mesh, parameters, sol);
    
%% End of simulation

    print_status(parameters, sol, 9);
