function C = project(a,b)
%PROJECT Calculate projection of vector(s) A onto vector B
%
% A can be an N-by-2 or N-by-3 matrix of N vecors in R2 or R3, and B is the
% 1-by-2 or 1-by-3 vector onto which to project.
% 
% In non-vectorised form, the equation looks like this:
%   project = @(a,b) dot(a,b)*b/norm(b);  
%
% In vectorised form, NORM calculates the matrix norm, not the row-wise
% length, so we calculate the length of the vectors explicitly
%
% USAGE:
%   C = project(X,V)
%   
%
% EXAMPLES:
%   C = project( rand(10,3), rand(1,3) )

% Joshua Martin, 16-Feb-2012

assert(size(b,1)==1,'B must be 1-by-2 or 1-by-3')

A = a;
B = b(ones(size(A,1),1),:); % SAME AS: repmat(b, [size(A,1),1]);
lenC = dot(A, B, 2) ./ sqrt( sum(B .* B, 2) );  
C = bsxfun(@times, lenC, B);
