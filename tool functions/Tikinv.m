%Tikinv    computes the regularized inverse of a matrix.
%
% [ invA, inv_S2 ] = Tikinv( A, p )
%
% outputs:
% invA is a n by m matrix, the regularized inverse of A
% inv_s2 is a sequence of the squared regularized sigular values
% compensated_s is a numerically homogenized singular values
% 
% inputs:
% A is a m by n matrix
% p is the regularization parameter, default is 0.05
% 
% 
% 2017-2019 Szu-Yu Lee
% Bouma Lab - The Wellman Center for Photomedicine

function [ invA, inv_s2, compensated_s ] = Tikinv( A, p, varargin )
% A the original matrix to be inverted
% p controls the Tikhonov parameter
A_p = A';

[~, S2, V] = svd(A_p*A);                                                    
s = sqrt(diag(S2));
maxsin = max(s);

if nargin == 1
    p = 0.1;
end

Tiklambda = p*maxsin;

inv_s2 = 1./( s.^2 + Tiklambda^2 );
compensated_s = s.*sqrt(inv_s2);
  
invA = (V .* (inv_s2.')) * V' * A_p;                                          
% A' already contain one Sigma, so inv_S2 should cancel one and leave one inverse

end

