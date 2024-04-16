function [figure2] = plot_solution(mesh, parameters, sol)     
 
    print_status(parameters, sol, 7);

     % Define used method for printing in legend.
     if parameters.space_time == 1
        if parameters.stab == 0;     
            method = 'space-time/Galerkin';
        elseif parameters.stab == 1;
            method = 'space-time/artificial diffusion';
        elseif parameters.stab == 2; 
            method = 'space-time/streamline-upwind';
        elseif parameters.stab == 3;
            method = 'space-time/SUPG';
        elseif parameters.stab == 4;
            method = 'space-time/GLS';
        else
            method = 'Unknown formulation';    
        end
     else       
        if parameters.stab == 0;     
            method = 'standard Galerkin';
        elseif parameters.stab == 1;
            method = 'artificial diffusion';
        elseif parameters.stab == 2; 
            method = 'streamline-upwind';
        elseif parameters.stab == 3;
            method = 'SUPG';
        elseif parameters.stab == 4;
            method = 'GLS';
        else
            method = 'Unknown formulation';      
        end
     end
    
if mesh.nsd == 1   
    % Initialize figure
    figure2 = figure;
    axes1 = axes('Parent',figure2);
    hold(axes1,'on');
    box(axes1,'on');
    
    % Define axis bounds depending on problem
    if strcmp(parameters.exact,'fig55')
        axis([0,1,-.1,1.5]);
    elseif strcmp(parameters.exact,'fig514')
        axis([0,1,-.1,1.2]);
    elseif strcmp(parameters.exact,'fig516')
        axis([0,1,-.2,0.8]);
    end
    
    % Print strings
    str = sprintf('FE-Solution: $$a =%0.3g,\\; \\nu=%.3g,\\; Pe= %0.3g,\\;C=%.3g$$',...
        parameters.advection(1), ...
        parameters.diffusion, ...
        parameters.Pe_x, ...
        parameters.C);
    str2 = sprintf('Physical time: %0.3g s',sol.phys_time); 


    
    % Print title legend and axis labels
    title(str,'FontWeight','bold','FontSize',14,'Interpreter','latex');      
    xlabel('$x$','FontSize',14,'Interpreter','latex');
    ylabel('$u$','FontSize',18,'Interpreter','latex');
    

    % Initialize timer 
    time = 0.0;
    
    %Print initial solution
    h_init = plot(mesh.coord,sol.U(:,1),'-.','LineWidth',1.5);
    
    % Print initial exact solution        
    h_exact = plot(sol.x_exact,sol.U_exact(:,1),'-.b','LineWidth',1.5);    

    % Print initial FE solution
    if parameters.space_time == 1
        h_curr = plot(mesh.coord((mesh.nx_nodes+1):end,1),...
                sol.U((mesh.nx_nodes+1):end,1),'-r','LineWidth',1.5);
    else    
        h_curr = plot(mesh.coord,sol.U(:,1),'-r','LineWidth',1.5);
    end
    
    
    % Print time
    if strcmp(parameters.u_init,'fig516');
        str2 = sprintf('Physical time: %1.3g  s',time); 
        h_text = annotation(figure2,'textbox',[0.15 0.87 0.22 0.04],...
                'String',{str2},'FontSize',10,'FitBoxToText','off');        
    else        
        str2 = sprintf('Physical time: %1.3g  s',time); 
        h_text = annotation(figure2,'textbox',[0.15 0.75 0.2 0.04],...
                'String',{str2},'FontSize',10,'FitBoxToText','off'); 
    end

    for i = 1:parameters.tot_t_iter

        % Set frames per second for animation     
        pause(1/parameters.fps)

        time = parameters.dt * i ;
        % Remove previous time solution
        delete(h_exact);
        delete(h_curr);
        delete(h_text);

        % Print exact solution        
        h_exact = plot(sol.x_exact,sol.U_exact(:,i+1),...
                '-.b','LineWidth',1.5);    

        % Print FE solution
        if parameters.space_time == 1
            h_curr = plot(mesh.coord((mesh.nx_nodes+1):end,1),...
                    sol.U((mesh.nx_nodes+1):end,i+1),'-r','LineWidth',1.5);
        else
            h_curr = plot(mesh.coord,sol.U(:,i+1),'-r','LineWidth',1.5);
        end 
        
        % Print legend and time
        if strcmp(parameters.u_init,'fig516');

            str2 = sprintf('Physical time: %1.3g  s',time); 
            h_text = annotation(figure2,'textbox',[0.15 0.87 0.22 0.04],...
                    'String',{str2},'FontSize',10,'FitBoxToText','off'); 

            legend('Initial solution','Exact Solution',method, ...
                'Location','southwest');   

        else        
            str2 = sprintf('Physical time: %1.3g  s',time); 
            h_text = annotation(figure2,'textbox',[0.15 0.75 0.18 0.04],...
                    'String',{str2},'FontSize',10,'FitBoxToText','off'); 

            legend('Initial solution','Exact Solution',method, ...
                'Location','northwest');
        end


    end

elseif mesh.nsd ==2
                 
        % Reshaping data for plotting
        x = reshape(mesh.coord(:,1),mesh.nx_nodes,mesh.ny_nodes)';
        y = reshape(mesh.coord(:,2),mesh.nx_nodes,mesh.ny_nodes)';
        u = reshape(sol.U(:,1),mesh.nx_nodes,mesh.ny_nodes)';

        % Title sting
        str = sprintf('$$\\vec{a} = (%0.3g,%0.3g)^T,\\; \\nu=%.3g,\\; Pe_{norm} = %0.4g$$',...
                parameters.advection(1), ...
                parameters.advection(2), ...
                parameters.diffusion, ...
                parameters.Pe_norm);
            
        figure2 = figure;
        colormap('jet');
        axes1 = axes('Parent',figure2,'xTick',0:0.2:1.0,'yTick',0:0.2:1.0);
        hold(axes1,'on');
        view(axes1,[110,20]);
        grid(axes1,'on');
        h_curr = surf(x,y,u,'Parent',axes1);
        xlabel('$x$','FontSize',14,'Interpreter','latex');
        ylabel('$y$','FontSize',14,'Interpreter','latex');
        zlabel('$u$','FontSize',18,'Interpreter','latex');

        
        title(['Formulation: ',method,', ',str],...
            'FontWeight','bold',...
            'FontSize',14,...
            'Interpreter','latex');
        colorbar('peer',axes1);
        
        time = 0.0;
        
        str2 = sprintf('Physical time: %1.3g  s',time); 
        h_text = annotation(figure2,'textbox',[0.15 0.75 0.18 0.04],...
                'String',{str2},'FontSize',10,'FitBoxToText','off'); 
            
        
        for i = 1:parameters.tot_t_iter

        u = reshape(sol.U(:,i),mesh.nx_nodes,mesh.ny_nodes)';
            
        % Set frames per second for animation     
        pause(1/parameters.fps)

        % Current time
        time = parameters.dt * i ;
        
        % Remove previous time solution
        delete(h_curr);
        delete(h_text);
 

        % Print FE solution
        h_curr = surf(x,y,u,'Parent',axes1);

        % Print legend and time      
        str2 = sprintf('Physical time: %1.3g  s',time); 
        h_text = annotation(figure2,'textbox',[0.15 0.75 0.18 0.04],...
                'String',{str2},'FontSize',10,'FitBoxToText','off'); 

        end


end

    print_status(parameters, sol, 8);
    
    
    
end
