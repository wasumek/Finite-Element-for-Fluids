function sol = exact_solution(mesh, parameters, sol)
% This function generates the exact solutions. Incase of the problems in
% Fig 5.5 and 5.14 of the FEF book, the solutions are read in from txt
% files present in the solver directory. 

if mesh.nsd == 1;
     
    % Read in solution from txt-files. 
    if ischar(parameters.exact) && (strcmp(parameters.exact,'fig55') ||...
            strcmp(parameters.exact,'fig514') )
    
        tmax = parameters.tot_t_iter * parameters.dt;
        filename = sprintf('exact_%.3gs_pe%.3g.txt',tmax,parameters.Pe_x);
        
        if exist(filename, 'file') == 2
            % Store read-in solution in sol.u_exact if it exists
            sol.u_exact = dlmread(filename);
            sol.x_exact = linspace(mesh.x0,mesh.x1,length(sol.u_exact));
            
        else
            % Set solution in sol.u_exact to zero if it doesn't exist
            sol.u_exact = zeros(10,1);    
            sol.x_exact = linspace(mesh.x0,mesh.x1,10)';
        end 
        
    elseif ischar(parameters.exact) && strcmp(parameters.exact,'fig516')
        % Generate exact solution for Gaussian hill problem.
        sol.x_exact = linspace(mesh.x0,mesh.x1,200)';
        x0 = 2/15;
        l = 7 * sqrt(2) / 300 ;

        sigma = sqrt( 1 + 4 * parameters.diffusion * sol.phys_time / ...
            (l^2) );
        
        sol.u_exact = 5 / 7 / sigma * exp( - ( ( sol.x_exact - x0 - ...
            parameters.advection(1) * sol.phys_time ) / (l * sigma) ).^2  );

    else
        % Set exact solution to zero if no exact solution is given.
        sol.u_exact = zeros(mesh.n_nodes,1);
        sol.x_exact = linspace(mesh.x0,mesh.x1,mesh.n_nodes);
    end
  
elseif mesh.nsd == 2;
    
    % For 2D case no exact solution is available. Hence the 
    % function's output are dummy variables.

    sol.u_exact = 0;
    sol.x_exact = 0;

end

end
