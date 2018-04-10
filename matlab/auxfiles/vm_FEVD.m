% Compute Mean Squared Errors of Forecasts up to period H;

function FEVD = vm_FEVD(P,phi,nlags_,H)

nv = size(P,1);
FEVD = zeros(nv,H);
variance = zeros(nv,H);   % Temporary array to store variances;

% Create Companion Form
F          = zeros(nv*nlags_,nv*nlags_);    % Matrix for Companion Form
I          = eye(nv);

for i=1:nlags_-1
    F(i*nv+1:(i+1)*nv,(i-1)*nv+1:i*nv) = I;
end
phi = phi(1:nv*nlags_,:);
F(1:nv,1:nv*nlags_)    = phi';


C = zeros(nv*nlags_,nv);
C(1:nv,:) = P;

for ii = 1:nv
    % Matrix that identifies shocks
    Ws = zeros(nv); 
    Ws(ii,ii)=1;
    v1 = C*Ws*C';
    vari = v1;
    
    % Store FEVD
    FEVD(:,1,ii) = diag(vari(1:nv,1:nv));
    
    for tt = 2:H
        vari = v1 + F*vari*F';
        FEVD(:,tt,ii) = diag(vari(1:nv,1:nv));
    end 
    variance = variance + squeeze(FEVD(:,:,ii));
end
       
for jj=1:nv        % For each variable
    for ii=1:nv    % For each shock
        FEVD(jj,:,ii)=100*squeeze(FEVD(jj,:,ii))./variance(jj,:);
    end
end
