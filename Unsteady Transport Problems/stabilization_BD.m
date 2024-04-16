function [Se_iq, Se_rhs_iq] = stabilization_BD(mesh,parameters,shape,integrator,ue_curr_t)
% This function assembles the stabilization element matrix for the 
% element under consisideration. 
%
% The stabilization method implemented here is the balancing diffusion
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

    
    if mesh.nsd == 1;
                
        beta = coth(parameters.Pe_x) - 1 / parameters.Pe_x;
        
        nu_bar = beta * parameters.advection(1) * mesh.h / 2; 

        Se_iq = parameters.theta * nu_bar * integrator * ...
                     (shape.Nx_iq' * shape.Nx_iq ...
                    + shape.Ny_iq' * shape.Ny_iq );

        Se_rhs_iq = - nu_bar * integrator * ...
                     (shape.Nx_iq' * shape.Nx_iq ...
                    + shape.Ny_iq' * shape.Ny_iq ) * ue_curr_t;        
    elseif mesh.nsd == 2;

        betaX = coth(parameters.Pe_x) - 1 / parameters.Pe_x;
        betaY = coth(parameters.Pe_y) - 1 / parameters.Pe_y; 

        nu_bar_ij = (parameters.advection(1) * betaX * mesh.hx ...
                   + parameters.advection(2) * betaY * mesh.hy)/2;

        nu_bar = nu_bar_ij / norm( parameters.advection )^2;
        
        Se_iq = parameters.theta * integrator * nu_bar * ( ...
          ( parameters.advection(1) * shape.Nx_iq ...
          + parameters.advection(2) * shape.Ny_iq)' ...
        * ( parameters.advection(1) * shape.Nx_iq ...
          + parameters.advection(2) * shape.Ny_iq) ) ;
      
        Se_rhs_iq = - integrator * nu_bar * ( ...
          ( parameters.advection(1) * shape.Nx_iq ...
          + parameters.advection(2) * shape.Ny_iq)' ...
        * ( parameters.advection(1) * shape.Nx_iq ...
          + parameters.advection(2) * shape.Ny_iq) ) * ue_curr_t;        

    end

   
end
