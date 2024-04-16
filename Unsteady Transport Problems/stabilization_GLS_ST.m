function [Se_iq, Se_rhs_iq] = stabilization_GLS_ST(mesh,parameters,shape,integrator, s_star)
% This function assembles the stabilization element matrix for the 
% element under consisideration. 
%
% The stabilization method implemented here is the Galerkin least squares
% method.
%
% Input:
% mesh              --> mesh data
% parameters        --> computational parameters
% shape             --> shape function data
% integrator        --> quadweight and Jacobian determinant for integration
%
% Output:
% Se_iq(n_elm_nodes,n_elm_nodes) -->     Fully assembled element matrix at
%                                        quadrature point iq.
% Se_rhs_iq(n_elm_nodes,1)       -->     Fully assembled element matrix at 
%                                        current quad point (RHS 
%                                        contribution).
%

        nu_bar = (( 1 / parameters.dt) ...
                 +( 2 * parameters.advection(1) / mesh.h ) ...
                 +( 4 * parameters.diffusion / mesh.h / mesh.h ) )^(-1);



        Pw = shape.Nt_iq' + ( parameters.advection(1) * shape.Nx_iq )';
        
        Ru1 = shape.Nt_iq + parameters.advection(1) * shape.Nx_iq;
           
        Se_iq = integrator * Pw * nu_bar * Ru1;
        
        Ru2 =  shape.N_iq * s_star ;
        
        Se_rhs_iq = integrator * Pw * nu_bar * Ru2;
        
end
   
