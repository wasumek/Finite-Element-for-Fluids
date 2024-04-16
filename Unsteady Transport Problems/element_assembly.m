function [ Ae, Fe ] = element_assembly( mesh, parameters, quads, shape,  x, ue_curr_t)

% This function assembles the element matrix of the element under
% consisideration. The contributions to an element level matrix
% (Ae(n_elm_nodes,n_elm_nodes)) are typically diffusive and advective terms
% following from the weak form. 
%
% Input:
% mesh                          --> Mesh data
% parameters                    --> Computational parameters
% quads                         --> Quadrature data
% shape                         --> Shape function data
% x(n_elm_nodes,nsd)            --> Nodal coordinates of current element.
% ue_curr_t (n_elm_nodes,nsd)   --> Nodal solution (u^n) of current element.
%
% Output:
% Ae(n_elm_nodes,n_elm_nodes) --> Fully assembled element matrix
% Fe(n_elm_nodes,n_elm_nodes) --> Fully assembled element rhs vector
%

%% Allocating the Diffusion (Ke) and Forcing or RHS (Fe) arrays.

    % LHS Contributions
	Te = zeros(mesh.n_elm_nodes,mesh.n_elm_nodes);      % Time der. term
    Ke = zeros(mesh.n_elm_nodes,mesh.n_elm_nodes);      % Diffusion
    Ce = zeros(mesh.n_elm_nodes,mesh.n_elm_nodes);      % Convection
    Se = zeros(mesh.n_elm_nodes,mesh.n_elm_nodes);      % Stabilization
    
    % RHS contributions
    Se_rhs = zeros(mesh.n_elm_nodes,1);                 % Stabilization
    Fe_rhs = zeros(mesh.n_elm_nodes,1);                 % RHS vector
    

%% Element level assembly
% Next the element matrix and rhs vector are assembled. This is achieved by
% a loop over the quadrature points in order to integrate over the
% elements.

    for iq = 1:quads.npoints;    

        % Compute shape function derivatives in physical space, the
        % Jacobians and their determinants.
        shape = shape_eval(mesh, parameters, shape, x, iq); 

        %% Setup integrator.
        % The integrator contains the quadrature weights and Jacobian
        % determinant which is responsible for the mapping. 

        integrator = quads.weights(iq) * shape.Jacobian_det;               

    	% Retrieve nodal values of current elements for the source term
    	[ s_star ] = gen_rhs( parameters, mesh);

        %% Assemble the element mass matrix 
        % This array contains the term dNa/dx * dNb/dx from the weak form. 
        
        Te = Te + shape.N_iq' * shape.N_iq / parameters.dt * integrator;
        
        %% Assemble the element stiffness matrix
        % This array contains the term dNa/dx * dNb/dx from the weak form.
        
        Ke = Ke + integrator * parameters.theta * parameters.diffusion *...
            (shape.Nx_iq' * shape.Nx_iq + shape.Ny_iq' * shape.Ny_iq );
                
        %% Assemble the element convection matrix Ce 
        % This array contains the term dNa * dNb/dx from the weak form.
        
    	Ce = Ce + integrator * parameters.theta * shape.N_iq' * ...
        	( parameters.advection(1) * shape.Nx_iq ...
        	+ parameters.advection(2) * shape.Ny_iq); 
                
        %% Assemble the element RHS forcing vector.
        Fe_rhs = Fe_rhs + shape.N_iq' *shape. N_iq * s_star * integrator;           

        %% Stabilization Matrix 
        if parameters.stab == 0;
        % No stabilization
            Se = Se;
            Se_rhs = Se_rhs;              
        elseif parameters.stab == 1;
        % Artificial diffusion        
            [Se_iq, Se_rhs_iq] = stabilization_AD(mesh,parameters,shape,...
                                integrator,ue_curr_t);
            Se = Se + Se_iq;
            Se_rhs = Se_rhs + Se_rhs_iq;
        elseif parameters.stab == 2;
        % Balancing diffusion
            [Se_iq, Se_rhs_iq] = stabilization_BD(mesh,parameters,shape,...
                                integrator,ue_curr_t);
            Se = Se + Se_iq;
            Se_rhs = Se_rhs + Se_rhs_iq;
        elseif parameters.stab == 3;
        % SUPG
            [Se_iq, Se_rhs_iq] = stabilization_SUPG(mesh,parameters,...
                                shape,integrator, ue_curr_t, s_star);
            Se = Se + Se_iq;
            Se_rhs = Se_rhs + Se_rhs_iq;
        elseif parameters.stab == 4;
        % GLS
            [Se_iq, Se_rhs_iq] = stabilization_GLS(mesh,parameters,...
                                shape,integrator, ue_curr_t, s_star);
            Se = Se + Se_iq;
            Se_rhs = Se_rhs + Se_rhs_iq;
        end 
    end  
 
    %% Assemble element matrix and vector
    
    Ae = Te + Ke + Ce + Se;
    
    Fe = Fe_rhs - 1 / parameters.theta * (Ke + Ce) * ue_curr_t + Se_rhs;
       
end
     
     
     
     
     
     
     
     
    
