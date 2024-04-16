function quads = quadrature()
% Quadrature points, and weights as well as the number of quadrature points
% are obtained through this function. The structure is as follows:
%
%   STRUCTURE:              DESCRIPTION                     DIMENSION:
%   quad __ npoints     --> Number of quadrature points     (1)
%        |_ points      --> Quadrature point coordinates    (npoints,nsd)
%        |_ weights     --> Quadrature weights              (npoints,1)
%



    % 2D Bilinear element
%     if mesh.elem_type == 1
% 
%             
%         % Number of quad-points
%         quads.npoints = 4;
% 
%         % Quad points
%         quads.points = 1/sqrt(3) * [-1 -1; 1 -1; 1  1; -1  1]; 
%                   
%         % Weights
%         quads.weights = [ 1; 1; 1; 1]; 

    % 2D Biquadratic  (Lagrangian quadrilateral)
    %elseif mesh.elem_type == 2 
        
        % Number of quad-points
        quads.npoints = 9;

        % Quad points
        quads.points = sqrt(3/5) * [-1 -1; 0 -1; 1 -1; 
                                    -1  0; 0  0; 1  0; 
                                    -1  1; 0  1; 1  1];
  
                  
        % Weights
        w1 = 5/9; w2 = 8/9; 
        quads.weights = [w1*w1; w2*w1; w1*w1;
                         w1*w2; w2*w2; w1*w2;
                         w1*w1; w2*w1; w1*w1];
        
        
   % end
        
    
end