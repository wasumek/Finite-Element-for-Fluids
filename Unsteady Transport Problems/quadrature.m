function quads = quadrature( mesh, parameters )
% Quadrature points, and weights as well as the number of quadrature points
% are obtained through this function. The structure is as follows:
%
%   STRUCTURE:              DESCRIPTION                     DIMENSION:
%   quad __ npoints     --> Number of quadrature points     (1)
%        |_ points      --> Quadrature point coordinates    (npoints,nsd)
%        |_ weights     --> Quadrature weights              (npoints,1)
%
% In case of the space-time formulation the space-only quadrature data  is 
% stored using the following arrays
%   quad __ npoints_s   --> Number of quadrature points     (1)
%        |_ points_s    --> Quadrature point coordinates    (npoints,nsd)
%        |_ weights_s   --> Quadrature weights              (npoints,1)
%
% Note that for the space thime solution in 1D the space-time mesh has
% two dimensions, i.e. space and time! The space-only mesh used for the
% jump term is 1D only. 


    % 1D Linear element
    if mesh.nsd == 1;
        
        if parameters.space_time == 1
            
            % Number of quad-points space only
            quads.npoints_s = 1;
             % Quad points space only
            quads.points_s = 0.0;
            % Weights space only
            quads.weights_s = 2.0;
            
            % Number of quad-points space-time
            quads.npoints = 4;

            % Quad points space-time
            quads.points = 1/sqrt(3) * [-1 -1; 1 -1; 1  1; -1  1]; 
                  
            % Weights space-time
            quads.weights = [ 1; 1; 1; 1]; 
            
        else
            % Single point Gauss quadrature
            if parameters.gauss_int  == 1

                % Number of quad-points
                quads.npoints = 1;

                % Quad points 
                quads.points = 0;

                % Weights
                quads.weights = 2;    

            % 2 Point Gauss quadrature
            elseif parameters.gauss_int  == 2

                % Number of quad-points
                quads.npoints = 2;

                % Quad points 
                quads.points = 1/sqrt(3) * [-1; 1];  

                % Weights
                quads.weights = [1; 1];   

            end
        end

    % 2D Bilinear element
    elseif mesh.nsd == 2;   
            
        % Number of quad-points
        quads.npoints = 4;

        % Quad points
        quads.points = 1/sqrt(3) * [-1 -1; 1 -1; 1  1; -1  1]; 
                  
        % Weights
        quads.weights = [ 1; 1; 1; 1]; 
        
    end
        
    
end