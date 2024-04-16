function [ A,  F ] = global_assembly( mesh, parameters, sol)
% This function is responsible for generating the global system matrix and
% rhs vector for the current time-step. This is done by a loop over the 
% elements. 
%
% Within the element loop each element matrix is generated using the
% "element_assembly()" function. The resulting element matrix is then 
% stored in the correct location of the global system matrix. The same is 
% done for the rhs vector. 
%
% Output of this routine is the fully assembled system matrix and rhs
% vector. Note that boundary conditions still need to be applied.
%
% Input:
% mesh       --> Mesh data
% parameters --> Computational parameters
% sol        --> Solution data
%
% Output:
% A          --> Fully assembled global system matrix
% F          --> Fully assembled global rhs vector
%

    %% Allocating and initializing arrays.
    A = zeros(mesh.n_nodes,mesh.n_nodes);
    F = zeros(mesh.n_nodes);


    x = zeros(mesh.n_elm_nodes,mesh.nsd);      % Elm. level coordinates
    ue_curr_t = zeros(mesh.n_elm_nodes,1);     % Elm. level solution at t=n


    if parameters.space_time == 1
        ue_curr_t_s = zeros(mesh.n_elm_nodes_s,1); % Elm. level solution at 
                                                   % t=n
        x_s = zeros(mesh.n_elm_nodes_s,mesh.nsd);  % Elm. level coordinates 
    end                                            % Of the space-only mesh 

    %% Retrieve the required Gauss quadrature data

    quads = quadrature( mesh, parameters );

    %% Evaluate shape functions and derivatives at quadrature points.

    % Remark: this function deals with the derivatives w.r.t. the reference 
    % coordinates.

    shape = shape_functions(mesh, parameters, quads);

    %% Loop over all elements
    for i = 1:mesh.n_elements
    
        % Retriev the space-time nodal 
        % coordinates and current 
        % time solution of the current 
        % element and store it in 
        % x_s(n_elm_nodes,nsd) and
        % ue_curr_t_s.
        if parameters.space_time == 1
            for j = 1:mesh.n_elm_nodes_s;
                x_s(j,:) = mesh.coord(mesh.conn(i,j),:);           
                ue_curr_t_s(j,:) = sol.u_t(mesh.conn(i,2+j),1);     
            end 
        end

        % Retriev the nodal coordinates and 
        % current time solution of the
        % current element and store it 
        % in x(n_elm_nodes,nsd) and
        % ue_curr_t_.
        for j = 1:mesh.n_elm_nodes;
            x(j,:) = mesh.coord(mesh.conn(i,j),:);
            ue_curr_t(j,:) = sol.u_t(mesh.conn(i,j),1);
        end

        % Generate element level matrix and 
        % rhs vector of current element
    
        if parameters.space_time == 1;
            % Space-time 
            [ Ae, Fe ] = element_assembly_ST( mesh, parameters, quads, shape, ...
                x, x_s, ue_curr_t_s);
        else
            % Semi-discrete

            [ Ae, Fe ] = element_assembly( mesh, parameters, quads, shape,  ...
                x,ue_curr_t);        
        end
    
        % Loop over all element nodes and store 
        % the element matrix and vector entries 
        % at the correct place in the global 
        % system matrix and rhs vector. 
        for j = 1:mesh.n_elm_nodes
            for k = 1:mesh.n_elm_nodes

                % Retrieve global matrix entry
                % location of current element
                % matrix entry.
                m = mesh.conn(i,j);           
                n = mesh.conn(i,k);

                % Store element matrix entries 
                % in global matrix.
                A(m,n) = A(m,n) + Ae(j,k);


            end       
            % Store element rhs vector 
            % entries in global rhs vector.
            F(m) = F(m) + Fe(j);    
        end
    end

    %% Store matrix A in an efficient (sparse) way.

    A = sparse(A);

end
