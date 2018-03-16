function rvs = igrand(s, nu)


x = randn(nu, 1);
rvs = sqrt(nu*s / sum(x.^2));

end