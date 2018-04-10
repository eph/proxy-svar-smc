function [ ln_iwish_pdf_val ] = pdf_ln_iwish( Psi, nu, SIGMA )
%PDF_LN_IWISH Log of the pdf of inverse wishart distribution
%   Psi: a square positive definite "scale" matrix
%   nu:  degrees of freedom, need nu > pdim - 1 where pdim = length(Psi)

pdim                =   size(Psi, 1); %dimension of square Psi matrix

ln_iwish_pdf_val    =    (nu/2) * log( det( Psi ) ) ...
                       - (nu * pdim/2) * log(2) ...
                       - ln_gamma_mv( pdim, nu/2 ) ...
                       - (nu + pdim + 1)/2 * log( det( SIGMA ) ) ...
                       -  0.5 * trace( Psi / SIGMA ) ;

end