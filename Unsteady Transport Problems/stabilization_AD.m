function [Se_iq, Se_rhs_iq] = stabilization_AD(mesh,parameters,shape,integrator, ue_curr_t)
% This function assembles the stabilization element matrix for the 
% element under consisideration. 
%
% The stabilization method implemented here is the artificial diffusion
% method.
%
% Input:
% mesh              --> mesh data
% parameters        --> computational parameters
% shape             --> shape function data
% integrator        --> quadweight and Jacobian determinant for integration
% ue_curr_t         --> current time solution at element
%
% Output:
% Se_iq(n_elm_nodes,n_elm_nodes) -->     Fully assembled element matrix at 
%                                        current quad point (LHS 
%                                        contribution).
% Se_rhs_iq(n_elm_nodes,n_elm_nodes) --> Fully assembled element matrix at 
%                                        current quad point (RHS 
%                                        contribution).
%

    % Artificial diffusion coefficient nu_bar. 
    
    if mesh.nsd == 1;
        
        nu_bar = parameters.advection(1) * mesh.h / 2; 

    elseif mesh.nsd == 2;
        
        nu_bar = norm(parameters.advection) * mesh.h / 2;  
    
    end

        Se_iq = parameters.theta * nu_bar * integrator * ...
            (shape.Nx_iq' * shape.Nx_iq + shape.Ny_iq' * shape.Ny_iq );
        Se_rhs_iq = - nu_bar * integrator * ...
            (shape.Nx_iq' * shape.Nx_iq + shape.Ny_iq' * shape.Ny_iq ) * ue_curr_t;

        
   
end
