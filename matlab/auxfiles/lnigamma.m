function y = lnigamma(x,a,b)
% LNPDFIG(X,A,B)
%	calculates log INVGAMMA(A,B) at X

% 03/03/2002
% Sungbae An
y = log(2) - gammaln(b/2) + (b/2)*log(b*a^2/2) - ( (b+1)/2 )*log(x^2) - b*a^2/(2*x^2);
