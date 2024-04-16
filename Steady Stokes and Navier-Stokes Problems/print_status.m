function print_status(sol, status)

    if status == 1
        
        display( ' ')    
        display(       '=============================')
        display(       '     Start Solution loop     ')
        display(       '=============================')   
        
    elseif status == 2
        
        str1 = sprintf(' Increment: %.3g L_inf: %.3g',sol.counter,sol.convergence);                    
        disp(str1)
        
    elseif status == 3
        
        display(       '=============================')
        display(       '   Solution loop finished!   ')
        display(       '=============================')
    
    end
        
end