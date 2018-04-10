data {
  // settings for the VAR
  int n;			
  int T;                   
  int p;			// number of lags
  int cons;			// include a constant

  // data 
  matrix[T,n] y;
  matrix[T,n*p+1] x;

  // prior distribution hyperparameters
  cov_matrix[n] S;
  int df;
  matrix[n*p+cons,n*p+cons] XtXinv;    // (X'X)^{-1}
  vector[n*(n*p+cons)] phihat;

  // proxy svar
  vector[T] mm;
  real<lower=0.0> sigmanu;
}

parameters {
  cov_matrix[n] sigma;
  vector[n*(n*p+cons)] phivec;
  real beta;
  vector[n*n] norms;
}

model {
  int m = n*p + cons;
  int i1; 
  int j1;
  matrix[m*n,m*n] sigma_kron_XtXinv;
  matrix[(n*p+cons),n] phi;
  matrix[T,n] yhat;
  
  matrix[n,n] Q;
  matrix[n,n] A0;
  matrix[(n*p+cons),n] Ap;
  matrix[T,n] eps;
  matrix[n,n] sigma_cholesky;
  // p(Phi,Sigma) = p(Phi|Sigma)p(Sigma)
  sigma ~ inv_wishart(df, S);

  // Sigma kron (X'X)^{-1}
  for (j in 1:n) {
    for (i in 1:n) {
      i1 = (i-1)*(m) + 1;
      j1 = (j-1)*(m) + 1;
      sigma_kron_XtXinv[i1:i1+m-1,j1:j1+m-1] = sigma[i,j]*XtXinv;
    }
  }
  phivec ~ multi_normal(phihat, 0.5*sigma_kron_XtXinv+0.5*sigma_kron_XtXinv');

  // p(Y|Phi,Sigma)
  phi = to_matrix(phivec, m, n);
  //yhat = x*phi;
  //for (i in 1:T) y[i] ~ multi_normal(yhat[i],sigma);

  // p(Omega)
  norms ~ normal(0,1);
  Q = qr_Q(to_matrix(norms,n,n));
  sigma_cholesky = cholesky_decompose(sigma);

  A0 = (Q' / sigma_cholesky);
  Ap = phi * A0;

  eps = y*A0 - x*Ap;

  beta ~ normal(0,1);
  mm ~ normal(beta*eps[,1],sigmanu);
}

