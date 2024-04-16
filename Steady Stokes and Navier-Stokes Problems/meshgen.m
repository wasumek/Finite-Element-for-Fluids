function mesh = meshgen(x0, x1, y0, y1, nx_elements, ny_elements, nsd, parameters)
% This routine calls the 2D mesh generator for Q1Q1 and Q2Q1 elements
% As output the routine stores all the input data
% and mesh data in the structure "mesh"

    if parameters.elem_type == 1
        % Q1Q1 element
        mesh = mesh_2D_Q1Q1(x0, x1, y0, y1, nx_elements, ny_elements, nsd, parameters);
    elseif parameters.elem_type == 2 
        % Q2Q1 element
        [ mesh ] = mesh_2D_Q2Q1(x0, x1, y0, y1, nx_elements, ny_elements, nsd, parameters);
    end
    
end

