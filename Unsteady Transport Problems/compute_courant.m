function parameters = compute_courant(mesh, parameters)
% This function computes the Courant number for 1D or 2D. 
  
    % Compute Courant number
    if mesh.nsd == 1;
        parameters.C = parameters.advection(1) * parameters.dt / mesh.h;
    elseif mesh.nsd == 2;
        parameters.C = norm(parameters.advection) * parameters.dt / mesh.h;
    end
    
end
