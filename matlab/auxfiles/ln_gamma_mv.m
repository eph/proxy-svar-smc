function lnval = ln_gamma_mv(p, a)
% GAMMA_MV_LN log of multivariate gamma function
    
lnval = p*(p-1)/4 * log(pi) + sum( gammaln( a + (1 - (1:p)) / 2) );
end

