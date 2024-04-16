function elmat = element_assembly( mesh, parameters, quads, shapeP, shapeV, xP, xV, ue, ve )
% This function assembles the element matrix for the element under
% consisideration. The contributions to an element level matrix
% (Ae(n_elm_nodes,n_elm_nodes)) are typically diffusive and advective terms
% as well as stabilization terms following from the weak form. 
%
% Input:
% mesh                  --> mesh data
% parameters            --> computational parameters
% quads                 --> quadrature data
% shapeP                --> Pressure shape functions and Jacobian data
% shapeV                --> Velocity shape functions and Jacobian data
% xP(n_p_elm_nodes,nsd) --> Nodal coordinates of current pressure element.
% xV(n_v_elm_nodes,nsd) --> Nodal coordinates of current velocity element.
% ue(n_v_elm_nodes,nsd) --> Nodal u-velocity of current velocity element.
% ve(n_v_elm_nodes,nsd) --> Nodal v-velocity of current velocity element.
%
% Output:
% elmat                 --> Contains all element matrices and vectors,
%                           e.g., Ke, Ce, Ge, Fe, He.
%

%% Allocating the element matrices and vectors.

    % LHS element matrices
    elmat.Ke = zeros(mesh.edfU,mesh.edfU); % Viscous matrix 
    elmat.Ce = zeros(mesh.edfU,mesh.edfU); % Convective matrix
    elmat.Ge = zeros(mesh.edfU,mesh.edfP); % Discrete gradient operator
    elmat.Pe = zeros(mesh.edfP,mesh.edfP); % Pres. matrix (typically empty)
    
    % RHS system vectors
    elmat.Fe = zeros(mesh.edfU,1);         % RHS F vector
    elmat.He = zeros(mesh.edfP,1);         % RHS H vector

    elmat.Se = zeros(mesh.edfU,mesh.edfU);
    elmat.SGe = zeros(mesh.edfU,mesh.edfP);
    elmat.SGTe = zeros(mesh.edfP,mesh.edfU);
    elmat.SPe  = zeros(mesh.edfP,mesh.edfP); 
    elmat.SHe  = zeros(mesh.edfP,1);
    
    % Matrix / vector indices at which to store u,v,p contributions
    U = 1:1:mesh.n_v_elm_nodes;
    V = mesh.n_v_elm_nodes + 1:1:2*mesh.n_v_elm_nodes;
    P  = 1:1:mesh.n_p_elm_nodes;

%% Element level assembly
% Next the element matrices and rhs vectors are assembled. This is achieved 
% by a loop over the quadrature points. By adding all quadrature point
% contributions the integrals are directly computed here.

    for iq = 1:quads.npoints 

        % Compute shape functions and derivatives in physical space, the
        % Jacobians and their determinants. Shape functions on the pressure
        % and velocity elements are stored in shapeP and shapeV
        % respecitvely.
        shapeP = shape_eval(mesh, shapeP, xP, iq);       
        shapeV = shape_eval(mesh, shapeV, xV, iq);
        
        % Compute the current velocity at current the quadrature point.
        uq = shapeV.N_iq * ue;
        vq = shapeV.N_iq * ve;
        
       	% Retrieve nodal source term values of current elements.
        rhs = gen_rhs( parameters, xV); 
        f_x = shapeV.N_iq * rhs(:,1);
        f_y = shapeV.N_iq * rhs(:,2);
        
        %% Setup integrator.
        % The integrator contains the quadrature weights and Jacobian
        % determinant which is responsible for the mapping. 

        integrator = quads.weights(iq) * shapeV.Jacobian_det;
        
        integrator_visc = parameters.visc * integrator; 

        %% Assemble the element viscosity matrix Ke 
        %  See also equation 6.24 of the book. 
        
        elmat.Ke(U,U) = elmat.Ke(U,U) + integrator_visc * ...
                ( shapeV.Nx_iq' * shapeV.Nx_iq ...
                + shapeV.Ny_iq' * shapeV.Ny_iq);        
            
        elmat.Ke(V,V) = elmat.Ke(V,V) + integrator_visc * ...
                ( shapeV.Nx_iq' * shapeV.Nx_iq ...
                + shapeV.Ny_iq' * shapeV.Ny_iq); 
            
        %% Assemble the element Convection matrix
        %  See also equation 6.33 of the book.
        
        if parameters.stokes ~= 1
            elmat.Ce(U,U) = elmat.Ce(U,U) + integrator * ...
                ( shapeV.N_iq' * shapeV.Nx_iq * uq ...
                + shapeV.N_iq' * shapeV.Ny_iq * vq );

            elmat.Ce(V,V) = elmat.Ce(V,V) + integrator * ...
                ( shapeV.N_iq' * shapeV.Nx_iq * uq ...
                + shapeV.N_iq' * shapeV.Ny_iq * vq );            
        end
        
        %% Assemble the element Pressure matrix G and G^t

        elmat.Ge(U,P) = elmat.Ge(U,P) - integrator * ( shapeV.Nx_iq' * shapeP.N_iq);
        elmat.Ge(V,P) = elmat.Ge(V,P) - integrator * ( shapeV.Ny_iq' * shapeP.N_iq);
        
        %% PSPG formulation
        % Only supported for Stokes problem.
        
        if (parameters.pspg == 1 && ...
            parameters.stokes == 1 && ...
            parameters.elem_type ==1 )
            
            tau = 1/3 * mesh.h_curr^2 / parameters.visc / 4;
       
            elmat.SPe(P,P) = elmat.SPe(P,P) - integrator * tau *...
                ( shapeP.Nx_iq' * shapeP.Nx_iq  + shapeP.Ny_iq' * shapeP.Ny_iq); 
            
            elmat.SHe(P) = elmat.SHe(P) - integrator * tau * ...
                ( shapeP.Nx_iq' * (shapeV.N_iq * rhs(:,1)) ...
                + shapeP.Ny_iq' * (shapeV.N_iq * rhs(:,2))); 
        end

        %% Assemble the element RHS vector Fe   
                
          	elmat.Fe(U) = elmat.Fe(U) + integrator * shapeV.N_iq' * f_x; 
            elmat.Fe(V) = elmat.Fe(V) + integrator * shapeV.N_iq' * f_y; 
            
    end  
    
end
     
     
     
     
     
     
     
     
    
