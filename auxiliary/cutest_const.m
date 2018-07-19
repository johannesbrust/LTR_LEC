function [varargout] = cutest_const( x, varargin )
  % 08/28/17
  % Assumption that CUTEst constraints are evaluated as: Ax - b0 = b
  
  if nargout > 3
      error( 'obj: too many output arguments' );
  end

  if nargin == 1
      % Compute objective function value
      if nargout > 1 
         % Gradient is requested
         
         [b,A]=cutest_cons(x);
         
         varargout{2} = b;
         varargout{1} = A;
        
         %[varargout{1}, varargout{2}] = cuter_obj(x);
      else
        varargout{1} =  cutest_cons(x);
      end
  else  
      % Only gradient is requested
      [b,A]             = cutest_cons(x);
      [varargout{1}]    = A;
  end