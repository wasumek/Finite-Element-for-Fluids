function sol = solve_time(mesh,parameters )
% This functions contains a loop over time which solved u_delta. This
% solution is added to the total solution as follows:
%
%   u_n+1 = u_n + u_delta
%
% The solutions u_n+1 are stored in the matrix sol.U, i.e. this matrix
% contains the time history. 



%% Initialize timestepping related arrays 

    % Current physical time
    sol.phys_time = 0;
    
    % Current time iteration
    sol.t_iter = 0;
    
    % Initialize current solution array with initial condition
    sol = initial_solution(mesh, parameters, sol);

    % Initialize first exact solution
    sol = exact_solution(mesh, parameters, sol); 
    sol.U_exact = sol.u_exact;
    
%% Time stepping loop

    % Print time step start
    print_status(parameters, sol, 1);

for i =1:parameters.tot_t_iter;
    
     % Update current physical time
    sol.phys_time = sol.phys_time + parameters.dt;
   
    % Update current iteration number
    sol.t_iter = i;   
    
    % Print current time-step data
    print_status(parameters, sol, 2);  
    
  %% Global Assembly
% Next the matrix and rhs vector of the linear system Au=f are generated.
% It is strongly suggested to study what happens inside this routine. The
% general steps in the routine in order to setup the matrix A and rhs 
% vector F are:
%
%   For each element:
%         * Setup the element matrix contribution.
%         * Setup the element rhs vector contribution.
%         * Store the element matrix and rhs vector correctly inside the 
%           global matrix.
%
%   The output of this routine is the global system matrix and vector "A" 
%   and "F" with size:
%
%   A(n_nodes,n_nodes)
%   F(n_nodes,1)
%

    % Assembling the global system matrices
    [ An, Fn ] = global_assembly( mesh, parameters, sol);


    %% Apply Boundary Conditions to An and Fn

    [An, Fn, sol] = apply_bc(mesh, parameters, sol, An, Fn, 1);

    %% Solve matrix problem
    sol.u_delta = An \ Fn(:,1);

    if parameters.space_time == 1
        
        % For space-time case we directly solve u_n+1 (instead of u_delta).
        % Hence no extra steps or BC need to be applied. 
        sol.u_t = sol.u_delta;
        
    else
        %% Apply boundary conditions to full time solution
        [An, Fn, sol] = apply_bc(mesh, parameters, sol, An, Fn, 2);

        %% Construct full time solution
        sol.u_t = sol.u_t + sol.u_delta;
        
    end

    % Store current time solution into solution history array.
    sol.U = [sol.U sol.u_t]; 

    %% Solve exact solution
    sol = exact_solution(mesh, parameters, sol);
    
    % Store exact current time solution  into exact solution history array.
    sol.U_exact = [sol.U_exact sol.u_exact];
    
end

    print_status(parameters, sol, 5);

% End of function
end
