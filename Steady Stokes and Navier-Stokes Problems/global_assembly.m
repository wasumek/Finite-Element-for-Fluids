function [ A,  F ] = global_assembly(mesh, parameters, bc, sol)
% This function is responsible for generating the global system matrix and
% rhs vector. This is done by a loop over the elements. 
%
% Within the element loop each element matrix is generated using the
% "element_assembly()" function. The resulting element matrix is then 
% stored in the correct location of the global system matrix. The same is 
% done for the rhs vector. 
%
% Input:
% mesh       --> Mesh data
% parameters --> Computational parameters
% bc         --> Boundary condition data
% sol        --> Solution data from previous iteration
%
% Output:
% A          --> Fully assembled global system matrix
% F          --> Fully assembled global rhs vector
%

    %% Allocating and initializing arrays.
    K =  spalloc(mesh.dfU,mesh.dfU,10*mesh.edfU); 
    G =  spalloc(mesh.dfU,mesh.dfP,10*mesh.edfU);
    GT = spalloc(mesh.dfP,mesh.dfU,10*mesh.edfU); % G^T
    P =  spalloc(mesh.dfP,mesh.dfP,10*mesh.edfU);
    F =  zeros(mesh.dfU,1);
    H =  zeros(mesh.dfP,1);

    
    %% Retrieve the required Gauss quadrature data

    quads = quadrature( );
    
    %% Evaluate shape functions and derivatives at quadrature points.

    % Remark: these functions deal with the derivatives w.r.t. reference 
    % coordinates.

    % Pressure related shape functions ( function space Q )
    shapeP = shape_functions(quads, mesh.p_elem_type);
    
    % Velocity related shape functions ( function space V/S )
    shapeV = shape_functions(quads, mesh.v_elem_type);
    

    %% Loop over all elements
    for i = 1:mesh.n_elements

        % Retrieve the nodal coordinates of 
        % the current element and store it 
        % temporarily in xP(n_p_elm_nodes,nsd)
        % and xV(n_v_elm_nodes,nsd) for pressure
        % and velocity respectively.
        xP = mesh.Pcoord(mesh.Pconn(i,:),:);
        xV = mesh.Vcoord(mesh.Vconn(i,:),:);
        
        % Retrieve the nodal velocity components 
        % of the current element and store it 
        % temporarily in ue(n_v_elm_nodes,nsd)
        % and ve(n_v_elm_nodes,nsd).
        M = 2*mesh.Vconn(i,:)-1;
        N = 2*mesh.Vconn(i,:);
        ue = sol.u(M);
        ve = sol.u(N);
        
        % Retrieve current mesh size
        mesh.h_curr = mesh.h_elm(i);

        % Generate element level matrix and 
        % rhs vector of current element
        elmat = element_assembly( mesh, parameters, quads, shapeP, shapeV, xP, xV, ue, ve );


        % Store the element matrices and vectores 
        % in their global equivalent. RS is used to 
        % define the global node numbering of u and v 
        % while PP is used to define the global node 
        % numbering of p.
        
        RS = [2*mesh.Vconn(i,:)-1 2*mesh.Vconn(i,:)];
        PP = mesh.Pconn(i,:);
        
        K(RS,RS)  = K(RS,RS)  + elmat.Ke  + elmat.Ce;
        G(RS,PP)  = G(RS,PP)  + elmat.Ge  + elmat.SGe; 
        GT(PP,RS) = GT(PP,RS) + elmat.Ge' + elmat.SGTe;
        P(PP,PP)  = P(PP,PP)  + elmat.SPe;
        F(RS)     = F(RS)     + elmat.Fe;  
        H(PP)     = H(PP)     + elmat.He  + elmat.SHe;
    end
    
    % Combine the K,C,G and P matrix to generate the full system matrix and
    % do the same for the rhs vector using F and H. The result should look
    % like:   
    %       _          _         _ _
    %       | K + C  G   |      | F |
    %   A = |            |  F = |   |
    %       |_G^T     0 _|      |_H_| 
    %
    % Note that addiitional matrices are added to account for the boundary
    % conditions (bc.Accd and bc.bccd).
  
    
  	A = [K        bc.Accd'                            	G;
         bc.Accd  zeros(bc.n_Dir_nodes,bc.n_Dir_nodes)	zeros(bc.n_Dir_nodes,mesh.dfP);
     	 G'       zeros(mesh.dfP,bc.n_Dir_nodes)        P];

    F = [F ;	bc.bccd ; H];

    %% Store matrix A in an efficient (sparse) way.
    A = sparse(A);
    
end
