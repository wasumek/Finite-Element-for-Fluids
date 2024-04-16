function [ Ae, Fe ] = element_assembly_ST( mesh, parameters, quads, shape,  x, x_s, ue_curr_t_s)
% This function assembles the element matricex for the element under
% consisideration. The contributions to an element level matrix
% (Ae(n_elm_nodes,n_elm_nodes)) are typically diffusive and advective terms
% following from the weak form. 
%
% Input:
% mesh               --> mesh data
% parameters         --> computational parameters
% x(n_elm_nodes,nsd) --> Nodal coordinates of current element.
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
    Je_s = zeros(mesh.n_elm_nodes_s,mesh.n_elm_nodes_s);% Jump term

    
    % RHS contributions
    Fe_rhs = zeros(mesh.n_elm_nodes,1);                 % RHS vector
    Se_rhs = zeros(mesh.n_elm_nodes,1);                 % RHS stab. vector
    Je_rhs_s = zeros(mesh.n_elm_nodes_s,1);             % Jump term 

%% Element level assembly
% Next the element matrix and rhs vector are assembled. This is achieved by
% a loop over the quadrature points in order to integrate over the
% elements.

    % This first loop is over the quadpoints in space only. This integrates
    % the jump term over the domain Omega (instead of the full space-time
    % domain Q).
    for iq = 1:quads.npoints_s
        %evaluate the shape functions, derivatives and Jacobian
       shape = shape_eval_space(mesh, parameters, shape, x, x_s, iq);
       
        %% Setup integrator.
        % The integrator contains the quadrature weights and Jacobian
        % determinant which is responsible for the mapping. 
        
       integrator = quads.weights_s(iq) * shape.Jacobian_det_s;
       
       %% Construct the LHS of the element Jump term.
       
       Je_s = Je_s + integrator*shape.N_iq_s' * shape.N_iq_s;
       
       %% Construct the RHS of the element Jump term.
       
       Je_rhs_s = Je_rhs_s + ...
           integrator*shape.N_iq_s' * shape.N_iq_s * ue_curr_t_s; 
       
    end

    % Loop over quadpoints of the space-time quad points. This will
    % integrate over the full space-time domain Q. 
    for iq = 1:quads.npoints;    

        % Compute shape function derivatives in physical space, the
        % Jacobians and their determinants.
        shape = shape_eval(mesh, parameters, shape, x, iq); 

        %% Setup integrator.
        % The integrator contains the quadrature weights and Jacobian
        % determinant which is responsible for the mapping. 

        integrator = quads.weights(iq) * shape.Jacobian_det;               

    	% Retrieve nodal values of current elements for the source term
    	s_star = gen_rhs( parameters, mesh);

        %% Assemble the element mass matrix 
        % This array contains the term dNa/dx * dNb/dx from the weak form.    
           
        Te = Te + shape.N_iq' * shape.Nt_iq * integrator;

        %% Assemble the element stiffness matrix 
        % This array contains the term dNa/dx * dNb/dx from the weak form.

        Ke = Ke +  integrator * parameters.diffusion * shape.Nx_iq' * shape.Nx_iq;        
           
        %% Assemble the element convection matrix Ce 
        % This array contains the term dNa * dNb/dx from the weak form.
    	Ce = Ce + integrator * shape.N_iq' * parameters.advection(1) * shape.Nx_iq;          

        %% Assemble the element RHS forcing vector.
        Fe_rhs = Fe_rhs + ...
            (shape.N_iq' * shape. N_iq ) * s_star * integrator; 

                %% Stabilization Matrix 
        if parameters.stab == 0;
        % No stabilization
            Se = Se;
            Se_rhs = Se_rhs;              
        elseif parameters.stab == 1;
        % Artificial diffusion      
            Se_iq = stabilization_AD_ST(mesh,parameters,shape,integrator);
            Se = Se + Se_iq;
        elseif parameters.stab == 2;
        % Balancing diffusion
            Se_iq = stabilization_BD_ST(mesh,parameters,shape,integrator);
            Se = Se + Se_iq;
        elseif parameters.stab == 3;
        % SUPG
            [Se_iq, Se_rhs_iq] = stabilization_SUPG_ST(mesh,parameters,...
                                shape,integrator, s_star);
            Se = Se + Se_iq;  
            Se_rhs = Se_rhs + Se_rhs_iq;
        elseif parameters.stab == 4;
        % GLS
            [Se_iq, Se_rhs_iq] = stabilization_GLS_ST(mesh,parameters,...
                                shape,integrator, s_star);
            Se = Se + Se_iq;
            Se_rhs = Se_rhs + Se_rhs_iq;
        end 
        
    end  
    %% Construct complete element level matrix Ae
    Ae = Te + Ke + Ce + Se;
    
    % Add jump term contributions to the lower time level only
    Ae(1:mesh.n_elm_nodes_s,1:mesh.n_elm_nodes_s) = ...
                Ae(1:mesh.n_elm_nodes_s,1:mesh.n_elm_nodes_s) + Je_s;

    %% Construct complete element level vector Fe
    Fe = Fe_rhs + Se_rhs;
    
    % Add jump term contributions to the lower time level only
    Fe(1:mesh.n_elm_nodes_s,1) = Fe(1:mesh.n_elm_nodes_s,1) + Je_rhs_s;
    
       
end
     
     
     
     
     
     
     
     
    
