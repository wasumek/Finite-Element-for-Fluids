function Se_iq = stabilization_AD_ST(mesh,parameters,shape,integrator)
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
%
% Output:
% Se_iq(n_elm_nodes,n_elm_nodes) --> Fully assembled element matrix at 
%                                    current quadpoint.
%

    % Artificial diffusion coefficient nu_bar. 

    nu_bar = parameters.advection(1) * mesh.h / 2; 

    Se_iq = nu_bar * integrator * shape.Nx_iq' * shape.Nx_iq;

        
   
end
