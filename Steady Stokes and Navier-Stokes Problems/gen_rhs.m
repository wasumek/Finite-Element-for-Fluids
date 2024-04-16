function [ rhs ] = gen_rhs( parameters, X )
% This function generates RHS source term at the current element
    
    x = X(:,1);
    y = X(:,2);

    aux_x = 1+(-4+(12-8.*y).*y).*y...
             +(-2+(24+(-72+48.*y).*y).*y ...
             +(12+(-48+(72-48.*y).*y).*y...
             +(-24+48.*y+(12-24.*y).*x).*x).*x).*x;
         
    aux_y = -12.*y.^2.*(1-y).^2 ...
            +(4+(-24+(48+(-48+24.*y).*y).*y).*y...
            +(-12+(72-72.*y).*y...
            +(8+(-48+48.*y).*y).*x).*x).*x;    

    
    rhs(:,1) = parameters.s(1) * aux_x;
    rhs(:,2) = parameters.s(2) * aux_y;
    
end
