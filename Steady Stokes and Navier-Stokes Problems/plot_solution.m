function plot_solution(mesh, bc, sol, parameters)
% within this function exact (if available) and FE solutions are plotted.

    % Reshaping data for plotting        
    x = reshape(mesh.Vcoord(:,1),mesh.nx_v_nodes,mesh.ny_v_nodes)';
    y = reshape(mesh.Vcoord(:,2),mesh.nx_v_nodes,mesh.ny_v_nodes)';
   
    xP = reshape(mesh.Pcoord(:,1),mesh.nx_p_nodes,mesh.ny_p_nodes)';
    yP = reshape(mesh.Pcoord(:,2),mesh.nx_p_nodes,mesh.ny_p_nodes)';
    
    uv = reshape(sol.u(1:mesh.n_v_nodes*mesh.nsd),2,mesh.n_v_nodes);        
    us = reshape(uv(1,:),mesh.nx_v_nodes,mesh.ny_v_nodes)';
    vs = reshape(uv(2,:),mesh.nx_v_nodes,mesh.ny_v_nodes)';

    ps = sol.u(mesh.dfU+bc.n_Dir_nodes+1:end);

    ps = reshape(ps,mesh.nx_p_nodes,mesh.ny_p_nodes)';
    ps = ps - ones(mesh.nx_p_nodes,mesh.ny_p_nodes)' * ...
            ( ps(1,1) + parameters.pshift );
    
    
    if parameters.matrix_fig == true
        figure2 = figure;
        axes1 = axes('Parent',figure2);
        hold(axes1,'on');
        box(axes1,'on');

        spy(sol.A)

        xlabel('$ndf \times nn$','FontSize',18,'Interpreter','latex');
        ylabel('$ndf \times nn$','FontSize',18,'Interpreter','latex');
        str = sprintf('Non-zero entries in global matrix');
        title(str,'FontWeight','bold','FontSize',18,'Interpreter','latex');
    end 

   if parameters.uvp_fig == true
        figure3 = figure;

        % Title sting
        str1 = 'velocity u';
        str2 = 'velocity v';
        str3 = 'pressure p';
                str4 = 'Normalized velocity vectors';
               str5 = 'Sreamlines';

        colormap('jet');
        axes1 = axes('Parent',figure3,'xTick',0:0.2:1.0,'yTick',0:0.2:1.0);
        axes2 = axes('Parent',figure3,'xTick',0:0.2:1.0,'yTick',0:0.2:1.0);
        axes3 = axes('Parent',figure3,'xTick',0:0.2:1.0,'yTick',0:0.2:1.0);
        hold(axes1,'on');
        hold(axes2,'on');
        hold(axes3,'on');
        view(axes1,[0.1,90]);
        grid(axes1,'on');
        view(axes1,[0.1,90]);
        grid(axes2,'on');
        view(axes1,[0.1,90]);
        grid(axes3,'on');
        
        subplot1 = subplot(2,3,1,'Parent',figure3);      
        hold(subplot1,'on');
        box on; axis equal;
        h_curr1 = contour(x,y,us);
        xlabel('$x$','FontSize',14,'Interpreter','latex');
        ylabel('$y$','FontSize',14,'Interpreter','latex');
        title([str1],...
            'FontWeight','bold',...
            'FontSize',14,...
            'Interpreter','latex');

        subplot2 = subplot(2,3,2,'Parent',figure3);
        hold(subplot2,'on');
        box on; axis equal;
        h_curr2 = contour(x,y,vs);
        xlabel('$x$','FontSize',14,'Interpreter','latex');
        ylabel('$y$','FontSize',14,'Interpreter','latex');
        title([str2],...
            'FontWeight','bold',...
            'FontSize',14,...
            'Interpreter','latex');

        subplot3 = subplot(2,3,3,'Parent',figure3);
        hold(subplot3,'on');
        view(subplot3,[-44.3 21.2]);
        h_curr3 = surf(xP,yP,ps);
        xlabel('$x$','FontSize',14,'Interpreter','latex');
        ylabel('$y$','FontSize',14,'Interpreter','latex');
        title([str3],...
            'FontWeight','bold',...
            'FontSize',14,...
            'Interpreter','latex');

        subplot4 = subplot(2,3,4,'Parent',figure3);
        hold(subplot4,'on');
        L = sqrt(us.^2 + vs.^2); 
        %L = 1;
        quiver(x,y,us./L,vs./L);
        box on; axis equal;
        axis([mesh.x0,mesh.y1,mesh.y0,mesh.y1]);
        xlabel('$x$','FontSize',14,'Interpreter','latex');
        ylabel('$y$','FontSize',14,'Interpreter','latex');
        title([str4],...
            'FontWeight','bold',...
            'FontSize',14,...
            'Interpreter','latex');
        
        subplot5 = subplot(2,3,5,'Parent',figure3);
        L = sqrt(us.^2 + vs.^2); 
        sx = [0.03 0.15  0.2 0.5 0.9];
        sy = [0.03 0.15  0.2 0.5 0.1];
        streamline(x,y,us,vs,sx,sy)
        axis equal; axis([0,1,0,1]); box on;
        xlabel('$x$','FontSize',14,'Interpreter','latex');
        ylabel('$y$','FontSize',14,'Interpreter','latex');
        title([str5],...
            'FontWeight','bold',...
            'FontSize',14,...
            'Interpreter','latex');
      
    subplot6 = subplot(2,3,6,'Parent',figure3);     
    axis square                             
    title('Mesh with boundary nodes')       

    % plot mesh 
    patch('Faces', mesh.Pconn, ...
          'Vertices', mesh.Pcoord, ...
          'MarkerSize',10,'Marker','diamond', ...
          'FaceColor','w');

    hold on

    for i = 1:mesh.n_v_nodes
        x = mesh.Vcoord(i,1);
        y = mesh.Vcoord(i,2);
       plot(x,y,'Marker','o','Color','black',...
             'MarkerFaceColor','black','MarkerSize',5);
    end
        % plot boundary
    for i = 1:mesh.n_Dir_nodes
            x = mesh.Vcoord(mesh.boundary(i,2),1);
            y = mesh.Vcoord(mesh.boundary(i,2),2);
            plot(x,y,'Marker','o',...
                 'MarkerFaceColor','r','MarkerSize',5);
        if mesh.boundary(i,1) == 1
            x = mesh.Vcoord(mesh.boundary(i,2),1);
            y = mesh.Vcoord(mesh.boundary(i,2),2);
            plot(x,y,'Marker','o','Color','red',...
                 'MarkerFaceColor','r','MarkerSize',5);
        elseif mesh.boundary(i,1) == 2
            x = mesh.Vcoord(mesh.boundary(i,2),1);
            y = mesh.Vcoord(mesh.boundary(i,2),2);
            plot(x,y,'Marker','o','Color','blue',...
                 'MarkerFaceColor','blue','MarkerSize',5);
        elseif mesh.boundary(i,1) == 3
            x = mesh.Vcoord(mesh.boundary(i,2),1);
            y = mesh.Vcoord(mesh.boundary(i,2),2);
            plot(x,y,'Marker','o','Color','green',...
                 'MarkerFaceColor','green','MarkerSize',5);
        elseif mesh.boundary(i,1) == 4
            x = mesh.Vcoord(mesh.boundary(i,2),1);
            y = mesh.Vcoord(mesh.boundary(i,2),2);
            plot(x,y,'Marker','o','Color','yellow',...
                 'MarkerFaceColor','yellow','MarkerSize',5);
        end
    end
           
        
        
    end
    

end
