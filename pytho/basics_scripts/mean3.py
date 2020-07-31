function [varargout] = mean3(varargin)
% MEAN : Calculate the mean value.
% For vectors: returns a value.
% For matrices: returns the mean value of each column.
    for ii = 1:nargin 
        %ii
       x=varargin{ii};
       [ m , n ] = size ( x );
       if m == 1
         m = n ;
       end

       varargout(ii)={sum(x)/m};
    end
