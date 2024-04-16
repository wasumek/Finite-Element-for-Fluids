function rhs = gen_rhs( parameters, mesh)
% This function generates a specified RHS


    if mesh.nsd == 1;

        rhs = ones(mesh.n_elm_nodes,1) * parameters.s; 

    elseif mesh.nsd == 2;  

        rhs = ones(mesh.n_elm_nodes,1) * parameters.s;

    end   
    
end
