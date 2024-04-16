function [Se_iq, Se_rhs_iq] = stabilization_SUPG(mesh,parameters,shape,integrator, ue_curr_t, s_star)
% This function assembles the stabilization element matrix for the 
% element under consisideration. 
%
% The stabilization method implemented here is the SUPG method.
%
% Input:
% mesh              --> mesh data
% parameters        --> computational parameters
% shape             --> shape function data
% integrator        --> quadweight and Jacobian determinant for integration
% ue_curr_t         --> current time solution at element nodes
% s_star            --> source term at element nodes
%
% Output:
% Se_iq(n_elm_nodes,n_elm_nodes) -->     Fully assembled element matrix at 
%                                        current quad point (LHS 
%                                        contribution).
% Se_rhs_iq(n_elm_nodes,1)       -->     Fully assembled element matrix at 
%                                        current quad point (RHS 
%                                        contribution).
%

    % Artificial diffusion coefficient nu_bar. 

    if mesh.nsd == 1;
                
        nu_bar = (( 1 / parameters.theta / parameters.dt) ...
                 +( 2 * parameters.advection(1) / mesh.h ) ...
                 +( 4 * parameters.diffusion / mesh.h / mesh.h ) )^(-1);

        
    elseif mesh.nsd == 2;

        nu_bar = (( 1 / parameters.theta / parameters.dt)^2 ...
                 +( 2 * norm(parameters.advection) / mesh.h )^2 ...
                 +( 4 * parameters.diffusion / mesh.h / mesh.h )^2 )^(-1/2);
    end

        Pw = ( parameters.advection(1) * shape.Nx_iq ...
             + parameters.advection(2) * shape.Ny_iq)';
        
        Ru1 = shape.N_iq / parameters.dt + parameters.theta * ...
                ( parameters.advection(1) * shape.Nx_iq ...
                + parameters.advection(2) * shape.Ny_iq);
           
        Se_iq = integrator * Pw * nu_bar * Ru1;
        
        Ru2 = - ( parameters.advection(1) * shape.Nx_iq ...
                + parameters.advection(2) * shape.Ny_iq) * ue_curr_t ...
                + shape.N_iq * s_star;
        
        Se_rhs_iq = integrator * Pw * nu_bar * Ru2;
        
   
end
