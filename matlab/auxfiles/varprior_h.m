function [ydu, xdu, breakss] = varprior_h(nv,nlags,nex,hyp,premom)
% nv: number of variables
% nlags: number of lags
% nex: number of exogenous variables including intercept
% hyp: vector of hyperparameters
% premom: pre-sample moments [mean; std]

tau = hyp(1,:);
dd = hyp(2,:);
ww = hyp(3,:);
lam = hyp(4,:);
mm = hyp(5,:);

% initializations
dsize = nex + (nlags+ww+1)*nv;
breakss = zeros(5,1);
ydu = zeros(dsize,nv);
xdu = zeros(dsize,nv*nlags+nex);

% dummies for the coefficients of the first lag
sig = diag(premom(:,2));
ydu(1:nv,:) = tau*sig;
xdu(1:nv,:) = [tau*sig zeros(nv,(nlags-1)*nv+nex)];
breakss(1) = nv;

% dummies for the coefficients of the remaining lags
if nlags > 1
    ydu(breakss(1)+1:nv*nlags,:) = zeros((nlags-1)*nv,nv);
    j = 1;
    while j <= nlags-1
        xdu(breakss(1)+(j-1)*nv+1:breakss(1)+j*nv,:) = [zeros(nv,j*nv) tau*sig*((j+1)^dd) zeros(nv,(nlags-1-j)*nv+nex)];
        j = j+1;
    end % while
    breakss(2) = breakss(1)+(nlags-1)*nv;
else
    breakss(2) = breakss(1);
end % if

% dummies for the covariance matrix of error terms
ydu(breakss(2)+1:breakss(2)+ww*nv,:) = kron(ones(ww,1),sig);
breakss(3) = breakss(2)+ww*nv;

% dummies for the coefficents of the constant term
lammean = lam*premom(:,1)';
ydu(breakss(3)+1,:) = lammean;
xdu(breakss(3)+1,:) = [kron(ones(1,nlags),lammean) lam];
breakss(4) = breakss(3)+1;

% dummies for the covariance matrix of coefficients of different lags
mumean = diag(mm*premom(:,1));
ydu(breakss(4)+1:breakss(4)+nv,:) = mumean;
xdu(breakss(4)+1:breakss(4)+nv,:) = [kron(ones(1,nlags),mumean) zeros(nv,nex)];
breakss(5) = breakss(4)+nv;
