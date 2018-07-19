function [ f,g ] = object_quad_arg( x, c, Q )
%OBJECT_QUAD_ARG computes the function value and gradient
% of a quadratic objective function with symmetrix Hessian Q.

%{

    06/21/18 J.B.

    f = x'*c + 0.5* x'*Q*x,     g = c + Q*x 

%}

f   = x'*c + 0.5*x'*Q*x;

if nargout < 2
    return;
else
    g = c + Q*x;
end


end

