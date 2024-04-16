function [A, F, sol] = apply_bc(mesh, parameters, sol, A, F, flag)
% This routines modifies the A matrix and RHS vector to account for the
% boundary conditions. 
%
% Input:
% mesh       --> mesh data
% parameters --> computational parameters
% A          --> Global system matrix
% F          --> Globel system RHS vector
% flag       --> flag used for different stages in applying the BCs

% The next step is to set the boundary conditions. For the nodes at which
% dirichlet boundary conditions are set the solution does not need to be
% computed (hence the solution is known). In order to incorporate this
% information the global system matrix and rhs vector are slightly
% modified.
%

%% First stage
% For the global matrix A:
% For each i_th Dirichlet node, the corresponding row is set to
% zero. Additionally the i_th diagonal entry of the matrix is set to the
% identity. 
%
% For the global rhs vector:
% For the i_th Dirichlet node, the i_th entry in the vector is replaced by
% a zero. 


    if flag == 1

        if parameters.space_time == 1
        % Space-time BC
            for i = 1:mesh.n_nodes

                % Left boundary
                if mesh.boundary(i) == 1 && parameters.w1 == 1;

                    A(i,:) = 0;
                    A(i,i) = 1;
                    F(i,1) = parameters.u1;

                % Right boundary
                elseif mesh.boundary(i) == 2 && parameters.w2 == 1;

                    A(i,:) = 0;
                    A(i,i) = 1;
                    F(i,1) = parameters.u2;

                end

            end
            
        else
            
        % Semi-discrete BC
          	for i = 1:mesh.n_nodes


                % Right boundary
                if mesh.boundary(i) == 4 && parameters.w4 == 1;

                    A(i,:) = 0;
                    A(i,i) = 1;
                    F(i,1) = 0;

                % Bottom boundary
                elseif mesh.boundary(i) == 1 && parameters.w1 == 1;

                    A(i,:) = 0;
                    A(i,i) = 1;
                    F(i,1) = 0;

                % Left boundary
                elseif mesh.boundary(i) == 2 && parameters.w2 == 1;

                    A(i,:) = 0;
                    A(i,i) = 1;
                    F(i,1) = 0;

                % Top boundary
                elseif mesh.boundary(i) == 3 && parameters.w3 == 1;

                    A(i,:) = 0;
                    A(i,i) = 1;
                    F(i,1) = 0;

                end

            end
        end     
 %% Second stage
% For the global matrix A:
% For each i_th Dirichlet node, the corresponding row is set to
% zero. Additionally the i_th diagonal entry of the matrix is set to the
% identity. 
%
% For the solution vector sol.u_t:
% For the i_th Dirichlet node, the i_th entry in the vector is replaced by
% the prescribed Dirichlet value. 
    elseif flag == 2

        for i = 1:mesh.n_nodes

            % Right boundary
            if mesh.boundary(i) == 4 && parameters.w4 == 1;

                A(i,:) = 0;
                A(i,i) = 1;
                if mesh.coord(i,2) <= 0.2;
                    sol.u_t(i) = 0;
                else
                    sol.u_t(i) = parameters.u4;
                end

            % Bottom boundary
            elseif mesh.boundary(i) == 1 && parameters.w1 == 1;

                A(i,:) = 0;
                A(i,i) = 1;
                sol.u_t(i) = parameters.u1;

            % Left boundary
            elseif mesh.boundary(i) == 2 && parameters.w2 == 1;

                A(i,:) = 0;
                A(i,i) = 1;
                sol.u_t(i) = parameters.u2;

            % Top boundary
            elseif mesh.boundary(i) == 3 && parameters.w3 == 1;

                A(i,:) = 0;
                A(i,i) = 1;
                sol.u_t(i) = parameters.u3;

            end

        end

    end
end