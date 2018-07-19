function [A,b] = const_quad_arg(x,A,b)
% CONST_QUAD_2 constraint for a quadratic programming problem

    b  = A*x - b;

end