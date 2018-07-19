function [A,b] = const_quad_2(x)
% CONST_QUAD_2 constraint for a quadratic programming problem
% 07/17/17, J.B

    global dataf;

    A = dataf.A;
    
    b1 = dataf.b1;
    
    b  = A*x - b1;

end