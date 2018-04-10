function y = impulseresponce(F,J,Omega,H)

nv = size(Omega,1);
nlags_ = size(F,1)/nv;
Omega1 = [Omega;zeros((nlags_-1)*nv,size(Omega,2))];
y = J'*(F^(H-1))*Omega1;



