function output = variancedecompositionFD(A,J,Ssigmau,P,n,t,i_transf)
%This function computes the variance decomposition of a VAR model

ntransf = length(i_transf);
MSE1Aug = zeros(n+ntransf ,n+ntransf ); % initialize MSEdiag
nshock = size(P,2);
ylevTemp = zeros(ntransf,nshock);
ylevU = zeros(ntransf ,n);
W=zeros(n+ntransf ,1);
evp=zeros(nshock,n,n);
evp_aug=zeros(nshock,n+ntransf ,n+ntransf );
output = zeros(t+1,n+ntransf ,nshock);
Pchol = chol(Ssigmau,'lower');

for i=0:t 
    Pphi = (J'*(A^i)*J);
    Ttheta= Pphi*P;
    TthetaChol= Pphi*Pchol;
    ylevU = TthetaChol(i_transf,:) + ylevU;
    TthetaCholAug= [Pphi*Pchol; ylevU];
    MSE1Aug = TthetaCholAug*TthetaCholAug'  + MSE1Aug;
    MSEdiag = diag(MSE1Aug);

    for j = 1:nshock
          ylevTemp(:,j) = Ttheta(i_transf,j)+ ylevTemp(:,j);
          evp(j,:,:) =  Ttheta(:,j)*Ttheta(:,j)' + squeeze(evp(j,:,:));
          Ttheta_aug = [Ttheta(:,j); ylevTemp(:,j)];
          evp_aug(j,:,:) =  Ttheta_aug*Ttheta_aug' + squeeze(evp_aug(j,:,:));
          W(:,1)= diag(squeeze(evp_aug(j,:,:)))./MSEdiag;
          output(i+1,:,j)=W(:,1)' ;
    end
end
