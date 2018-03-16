function l =  loglik_m_given_y(m, u, sigma_tr, omega, bet, signu)

scale_mat = omega'/sigma_tr;


z = m' - bet*scale_mat*u';

l = sum(log(mvnpdf(z', [], signu.^2)));

end

