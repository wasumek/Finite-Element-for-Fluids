function [mesh] = meshgen(x0, x1, y0, y1, nx_elements, ny_elements, nsd, parameters)
% This routine calls the 1D or 2D mesh generator based on the number of
% spatial dimensions nsd. As output the routine stores all the input data
% and mesh data in the structure "mesh"

    if nsd == 1 
        
        if parameters.space_time == 1
            % Space-time mesh generation
            mesh = mesh_1D_ST(x0, x1, nx_elements, parameters);
        else
            % Conventional mesh generation
            mesh = mesh_1D( x0, x1, nx_elements, nsd, parameters);           
        end
    
    elseif nsd == 2;
        
        if parameters.space_time == 1
            % Space-time mesh generation not supported for 2D
            error('WARNING: Space-time formulation not supported for 2D problems');
        else 
            % Conventional mesh generation
            mesh = mesh_2D(x0, x1, y0, y1, nx_elements, ny_elements, nsd, parameters);
        end
        
    end

end

