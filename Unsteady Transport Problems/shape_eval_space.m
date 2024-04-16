function [ shape ] = shape_eval_space(mesh, parameters, shape, x, x_s, iq)
% This funcion returns the shape function values evaluated at the
% current quadrature point, the shape function derivatives in
% physical space, the Jacobian and its determinant. 
%
% The returned arrays are added to the "shape" structure and are 
% constructed as follows:
%
%   N_iq_s    = [N1(xi,eta), N2(xi,eta), ... , Nn(xi,eta)]
%
%   Nxi_iq_s  = [dN1/dxi(xi,eta), dN2/dxi(xi,eta), ... dNn/dxi(xi,eta)] 
%   Neta_iq_s = [dN1/deta(xi,eta), dN2/deta(xi,eta), ... dNn/deta(xi,eta)]
%
%   Nx_iq_s   = [dN1/dx(xi,eta), dN2/dx(xi,eta), ... dNn/dx(xi,eta)]
%   Ny_iq_s   = [dN1/dy(xi,eta), dN2/dy(xi,eta), ... dNn/dy(xi,eta)]
%   
%   Jacobian_s     = [dx/dxi   dy/dxi;
%                    dx/deta  dy/deta]
%   Jacobian_det_s = |J|
%
%   Note that xi and eta in the arrays above are those of the current quad
%   point 'iq'.
%

    %%   Retrieve the shape function values at the current quadpoint "iq"

    shape.N_iq_s = shape.N_s(iq,:);

    %% Shape function derivatives at current quadrature point. 
    % Here the derivatives of the shape functions in the reference domain
    % are retrieved for the current quatdature point. The arrays Nxi_iq 
    % and Neta_iq have dimensions (1,n_elm_nodes). When solving 1D problems
    % the Neta_iq array is automatically set to zero. 
    
    shape.Nxi_iq_s = shape.Nxi_s(iq,:);           

    % 2D Jacobian.
    shape.Jacobian_s = shape.Nxi_iq_s * x_s(:,1);

    
    %% Compute the Jacobian determinant |J|.
    
    shape.Jacobian_det_s = det(shape.Jacobian_s); 
    
    %% Compute the shape function derivatives in physical coordinates.
    % Obtain dN/dx and dN/dy by applying the chain rule to the shape 
    % functions by solving dN/dx = dNa/dxi * dxi/dx + dNb/dxi * dxi/dx.
    %
    % Use Matlab's built in equation solver. To solve Ax = b for x you can
    % simply write:   x = A\b
    
    shape.Nx_iq_s = shape.Jacobian_det_s * shape.Nxi_iq_s;
         
end


    