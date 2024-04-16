function [bc, sol] = solve_increment(mesh, parameters)
% This functions contains a loop over the increments. This
% loop solves u = u^* + u_delta until the increment solution has 
% reached convergence or the maximum number of iteration is reached.
% u^* in this case is the result of the previous iteration

%% Boundary conditions

    [bc, mesh] =  BC_Cavity(mesh, parameters);
 
%% Initialize the zeroth solution

    % reshape initial guess 
    u0 = reshape( [1;0] * ones(1,mesh.n_v_nodes), 1, mesh.dfU)' ... 
                        * parameters.u_init(1);
    v0 = reshape([0;1]*ones(1,mesh.n_v_nodes),   1,mesh.dfU)' ...
                        * parameters.u_init(2);   
   
    % Store the initial guess in the solution vector
    sol.u = [u0 + v0; zeros(bc.n_Dir_nodes,1); zeros(mesh.n_p_nodes,1)];

    % Initialize counter and residual
    sol.counter = 0;
	sol.convergence = 1;
    
%% Start iteration loop
	print_status(sol, 1);
    while ( sol.convergence > 1e-6 && sol.counter < parameters.n_iter );


        % Assembe the global system matrix and vector
        %
        % Here A and F contain the matrix and vector includint boumdary 
        % conditions:
        %    _          _    _ _
        %   | K + C  G^T |  | F |
        %   |            |  |   |
        %   |_  G     0 _|  |_H_| 
        %
        
        [ A,  F ] = global_assembly(mesh, parameters, bc, sol);
        
        % Compute LHS vector using the solution
        % of the previous iteration.
        aux = F - A * sol.u;
       
        % Compute the new increment solution
        u_delta = A \ aux;
        
        % Increment to new solution
        sol.u = sol.u + u_delta;
        
        % Compute Velocity residual and print the result
        sol.convergence = norm(u_delta(1:mesh.dfU,1),inf);
        sol.counter = sol.counter + 1;
        print_status(sol, 2);

   end
   sol.A = A;
   print_status(sol, 3);

% End of function
end
