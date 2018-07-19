function [ f,g ] = object_quad_arg_L( phi, x, c, Q )
%OBJECT_QUAD_ARG computes the function value and gradient
% of a quadratic objective function with symmetrix Hessian Q. For large
% problems.

%{

    06/21/18 J.B.

    f = x'*c + 0.5* x'*Q*x,     g = c + Q*x 


    07/02/18, J.B., large version

    Hessian = phi. I + Q Q'.

%}

f   = x'*c + (phi*0.5)*(x'*x) + 0.5*(x'*Q)*(Q'*x);

if nargout < 2
    return;
else
    g = c + phi.*x + Q*(Q'*x);
end


end

