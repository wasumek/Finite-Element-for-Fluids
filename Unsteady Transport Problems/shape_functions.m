function shape = shape_functions(mesh, parameters, quads)
% This functions computes the shape functions and their derivatives at a
% given quadrature point.
% 
% INPUT required:
% mesh          mesh data
% parameters    compuation parameters
% quads         quadrature data
% 
% OUTPUT:
%   STRUCTURE:          DESCRIPTION                 DIMENSION:
%   shape __ N      --> Shp fnct values at xi/eta   (npoints,n_elm_nodes)
%         |_ Nxi    --> Shp fnct derivatives to xi  (npoints,n_elm_nodes)
%         |_ Neta   --> Shp fnct derivatives to eta (npoints,n_elm_nodes)
%
% In case of the space-time formulation the space-only shape function data
% is stored using the following arrays
%
%   shape __ N_s    --> Shp fnct val at xi/eta  (npoints_s,n_elm_nodes_s)
%         |_ Nxi_s  --> Shp fnct der to xi      (npoints_s,n_elm_nodes_s)
%         |_ Neta_s --> Shp fnct der to eta     (npoints_s,n_elm_nodes_s)
%
%   The shape function arrays are constructed such that each row contains 
%   all shape functions evaluated at a specific quadrature point. This 
%   yields the following for the shape functions N1(xi),...Nn(xi) 
%
%                   Shape functions --->
%    _                                                   _
%   |  N1(xi_1,eta_1) N2(xi_1,eta_1) . . . Nn(xi_1,eta_1) |   |
%   |  N1(xi_2,eta_2) N2(xi_2,eta_2) . . . Nn(xi_2,eta_2) |   | 
%   |     .              .           .        .           |   | Quad points
%   |     .              .             .      .           |   |
%   |     .              .               .    .           |   V
%   |_ N1(xi_n,eta_n) N2(xi_n,eta_n) . . . Nn(xi_n,eta_n)_| 

%% Generate shape functions.
    
    % 1D linear elements.      
    if mesh.nsd == 1;
        
        % Space-time shape functions
        if parameters.space_time == 1
          
            %% Space only
             % Quadrature points coordinates in reference domain.
            xi_s = quads.points_s(:,1); 

            % Shape function values at quadrature point
            shape.N_s = [(1-xi_s)/2 (1+xi_s)/2]; 

            % Shape function xi-derivative values at quadrature points
            shape.Nxi_s = [-1/2 1/2]; 

            % Shape function eta-derivative values at quadrature points 
            % Thes are set to zero for 1D case. 
            shape.Neta_s = [0 0];
            
            %% Space-time
            % Quadrature points coordinates in reference domain.
            xi = quads.points(:,1); eta = quads.points(:,2);      

            % Shape function values at quadrature point
            shape.N = [ (1-xi).*(1-eta)/4 ...
                        (1+xi).*(1-eta)/4 ...
                        (1+xi).*(1+eta)/4 ...
                        (1-xi).*(1+eta)/4 ]; 

            % Shape function xi-derivative values at quadrature points
            shape.Nxi = [ (eta-1)/4 ...
                          (1-eta)/4 ...
                          (1+eta)/4 ...
                         -(1+eta)/4 ]; 

            % Shape function eta-derivative values at quadrature points
            shape.Neta = [ (xi-1)/4 ...
                          -(1+xi)/4 ...  
                           (1+xi)/4 ... 
                           (1-xi)/4 ];   
                       
        % Non space-time shape functions               
        else
      
            if parameters.gauss_int  == 1

                % Quadrature points coordinates in reference domain.
                xi = quads.points(:,1); 

                % Shape function values at quadrature point
                shape.N = [(1-xi)/2 (1+xi)/2]; 

                % Shape function xi-derivative values at quadrature points
                shape.Nxi = [-1/2 1/2]; 

                % Shape function eta-derivative values at quadrature points 
                % Thes are set to zero for 1D case. 
                shape.Neta = [0 0]; 

            elseif parameters.gauss_int  == 2;

                % Quadrature points coordinates in reference domain.
                xi = quads.points(:,1); 

                % Shape function values at quadrature point
                shape.N = [(1-xi)/2 (1+xi)/2]; 

                % Shape function xi-derivative values at quadrature points
                shape.Nxi = [-1/2 1/2 ;
                             -1/2 1/2 ]; 

                % Shape function eta-derivative values at quadrature points 
                % Thes are set to zero for 1D case. 
                shape.Neta = [0 0;
                              0 0]; 

            end

        end
        
    % 2D bilinear elements.
    elseif mesh.nsd == 2;
 
        % Quadrature points coordinates in reference domain.
        xi = quads.points(:,1); eta = quads.points(:,2); 
              
        % Shape function values at quadrature point
        shape.N = [ (1-xi).*(1-eta)/4 ...
                    (1+xi).*(1-eta)/4 ...
                    (1+xi).*(1+eta)/4 ...
                    (1-xi).*(1+eta)/4 ]; 

        % Shape function xi-derivative values at quadrature points
        shape.Nxi = [ (eta-1)/4 ...
                      (1-eta)/4 ...
                      (1+eta)/4 ...
                     -(1+eta)/4 ]; 
        
        % Shape function eta-derivative values at quadrature points
        shape.Neta = [ (xi-1)/4 ...
                      -(1+xi)/4 ...  
                       (1+xi)/4 ... 
                       (1-xi)/4 ]; 

    end
end