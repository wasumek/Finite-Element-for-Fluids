function print_status(parameters, sol, status);

    if status == 1;

        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        display(    '~ Start Time Stepping Loop! ~')
        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        display( ' ')
        display( ' ')
        tic;
    elseif status ==2;

        str1 = sprintf(' Time step    : %.3g of %.3g ',...
                            sol.t_iter,parameters.tot_t_iter);                    
        str2 = sprintf(' Physical time: %.3g         ',...
                            sol.phys_time);

        display( ' ')    
        display(       '=============================')
        disp(str1)
        disp(str2)
        display(       '=============================')
        display( ' ')   
    elseif status == 3

        str1 = sprintf(' Incr: %3.3g  L_inf: %8.3e  ',sol.n_iter,sol.residual);
        disp(str1)

    elseif status == 4

        str1 = sprintf(' Total Iterations: %.3g  ',sol.n_iter);
        str2 = sprintf(' L_inf residual:   %.3g  ',sol.residual);
        display(    '== Finished Inctement Loop ==')
        disp(str1)
        disp(str2)
        display(    '=============================')
        display( ' ') 

    elseif status == 5

        str1 = sprintf(' Time step: %4.3g of %4.3g ', sol.t_iter,parameters.tot_t_iter);
        display(    '== Finished Time Step Loop ==')
        disp(str1)
        display(    '=============================')
        display( ' ') 
        display( ' ') 

    elseif status == 6        
        time = toc;
        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        display(    '~~~ Computation Finished! ~~~')
        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp(sprintf(' Peclet number       :  %5.3g',parameters.Pe_x))
        disp(sprintf(' Peclet number y     :  %5.3g',parameters.Pe_y))
        disp(sprintf(' CFL number          :  %5.3g',parameters.C))
        disp(sprintf(' Space-time          :  %5.3g',parameters.space_time))
        disp(sprintf(' Time setp size      :  %5.3g',parameters.dt))
        disp(sprintf(' Computational time  : %4.3g',time))
        disp(' ')
        disp(' ')

    elseif status == 7       

        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        display(    '~~~ Postprocessing data! ~~~~')
        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp(' ')
        disp(' ')   
        
    elseif status == 8       

        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        display(    '~~ Postprocessing finished! ~')
        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp(' ')
        disp(' ')
    elseif status == 9       

        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        display(    '~~~~ FE Simulation done! ~~~~')
        display(    '~~~~         -           ~~~~')
        display(    '~~~~ Terminating solver. ~~~~')
        display(    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp(' ')
        disp(' ')

    end
end