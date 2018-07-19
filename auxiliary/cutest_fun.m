function [varargout] = cutest_fun( x, varargin )
  % Evaluate objective function
  % Usage:       f = obj(x)      evaluates function value only
  %          [f,g] = obj(x)  evaluates function value and gradient
  %        [f,g,H] = obj(x)  evaluates function value, gradient and Hessian

  if nargout > 3
      error( 'obj: too many output arguments' );
  end

  if nargin == 1
      % Compute objective function value
      if nargout > 1 
         % Gradient is requested
        [varargout{1}, varargout{2}] = cutest_obj(x);
      else
        varargout{1} =  cutest_obj(x);
      end
  else  
      % Only gradient is requested
      [varargout{1}] = cutest_grad(x);
  end