function parameters = compute_peclet(mesh, parameters)
% This function computes the Peclet number in 1D and 2D.
   
    % Compute Peclet number 
    parameters.Pe_x = parameters.advection(1) * mesh.hx / ...
                        (2 * parameters.diffusion);
                
    % Compute Peclet number 
    parameters.Pe_y = parameters.advection(2) * mesh.hy / ...
                        (2 * parameters.diffusion); 
    % Peclet_normed
    parameters.Pe_norm = norm(parameters.advection) * ...
                        max(mesh.hx, mesh.hy) / ...
                        (2 * parameters.diffusion);
   
    
end
