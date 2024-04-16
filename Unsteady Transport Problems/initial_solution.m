function sol = initial_solution(mesh,parameters, sol)
% Initialize FE and exact solution arrays sol.U and sol.U_exact.


    if parameters.u_init == 0

    sol.u_t = zeros(mesh.n_nodes,1);

    elseif strcmp(parameters.u_init,'fig516') 

        x0 = 2/15;
        l = 7*sqrt(2) / 300;
        
        sol.u_t = 5 / 7 * exp( - ( ( mesh.coord(:,1) - x0 ) / l ).^2  );
    end
    
    sol.U = sol.u_t;
    sol.U_exact = sol.u_t;
end