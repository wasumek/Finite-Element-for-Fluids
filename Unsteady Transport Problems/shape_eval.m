function [ shape ] = shape_eval(mesh, parameters, shape, x, iq)
% This funcion returns the shape function values evaluated at the
% current quadrature point, the shape function derivatives in
% physical space, the Jacobian and its determinant. 
%
% The returned arrays are added to the "shape" structure and are 
% constructed as follows:
%
%   N_iq    = [ N1(xi,eta), N2(xi,eta), ... , Nn(xi,eta)]
%
%   Nxi_iq  = [ dN1/dxi(xi,eta), dN2/dxi(xi,eta), ... , dNn/dxi(xi,eta)] 
%   Neta_iq = [ dN1/deta(xi,eta), dN2/deta(xi,eta), ... , dNn/deta(xi,eta)]
%
%   Nx_iq   = [ dN1/dx(xi,eta), dN2/dx(xi,eta), ... , dNn/dx(xi,eta)]
%   Ny_iq   = [ dN1/dy(xi,eta), dN2/dy(xi,eta), ... , dNn/dy(xi,eta)]
%   
%   Jacobian     = [dx/dxi   dy/dxi;
%                   dx/deta  dy/deta]
%   Jacobian_det = |J|
%
%   Note that xi and eta in the arrays above are those of the current quad
%   point 'iq'.
%
% When the space-time formulation is used similar arrays are added for the
% space-only case.
%

    %%   Retrieve the shape function values at the current quadpoint "iq"
    shape.N_iq = shape.N(iq,:);

    %% Shape function derivatives at current quadrature point. 
    % Here the derivatives of the shape functions in the reference domain
    % are retrieved for the current quatdature point. The arrays Nxi_iq 
    % and Neta_iq have dimensions (1,n_elm_nodes). When solving 1D problems
    % the Neta_iq array is automatically set to zero. 
        
    shape.Nxi_iq = shape.Nxi(iq,:);  
    shape.Neta_iq = shape.Neta(iq,:); 

    %% Compute the Jacobian J 
    % The Jacobian matrix for the 1D and 2D case are as follows
    %
    %   1D case:  J = dx/dxi
    %
    %                  
    %   2D case:  J = [dx/dxi  dy/dxi ;
    %                  dx/deta dy/deta]
    %
    % x and y have to be described as a function of the reference
    % coordinates. This yields for x and y in 2D:
    %
    %   x(xi,eta)  = sum ( N_i * x_i ),  
    %   y(xi,eta) =  sum ( N_i * y_i ), 
    % 
    % in this case only the shape functions are functions of the 
    % reference coordinates. Taking e.g. the derivative of
    % x(xi,eta) w.r.t. xi in 2D would yield:
    %
    % dx(xi,eta)/dxi = dN1/dxi * x1 + dN2/dxi * x2 + ... + dN4/dxi * x4   
    %
    % To compute the sum in the above example in an efficient way in Matlab
    % use can be made of simple matrix mulitiplication:
    %                                                    _  _
    %                                                   | b1 |
    %   a1*b1 + a2*b2 + ... + an*bn = [ a1 a2 a3 a4 ] * | b2 |
    %                                                   | b3 |
    %                                                   |_b4_|         
    
    if mesh.nsd == 1;
        % 1D Jacobian
        
        if parameters.space_time == 1
            % Space-time Jacobian
            shape.Jacobian = [shape.Nxi_iq * x(:,1) 0.0;   
                              shape.Neta_iq * x(:,1) 0.5 * parameters.dt];
        else 
            % Semi-discrete Jacobian.
            shape.Jacobian = shape.Nxi_iq * x(:,1);  
        end
        
    elseif mesh.nsd == 2;
        % 2D Jacobian.
        
        shape.Jacobian = [shape.Nxi_iq  * x(:,1) shape.Nxi_iq  * x(:,2);   
                          shape.Neta_iq * x(:,1) shape.Neta_iq * x(:,2)];
                      
    end
    
    %% Compute the Jacobian determinant |J|.
    
    shape.Jacobian_det = det(shape.Jacobian); 
    
    %% Compute the shape function derivatives in physical coordinates.
    % Obtain dN/dx and dN/dy by applying the chain rule to the shape 
    % functions by solving dN/dx = dNa/dxi * dxi/dx + dNb/dxi * dxi/dx.
    %
    % Use Matlab's built in equation solver. To solve Ax = b for x you can
    % simply write:   x = A\b
    
    tmp = shape.Jacobian\[shape.Nxi_iq;shape.Neta_iq]; 

    shape.Nx_iq = tmp(1,:); 
    
    if parameters.space_time == 1
        shape.Nt_iq = tmp(2,:);
    else
        shape.Ny_iq = tmp(2,:);   
    end
   
end


    