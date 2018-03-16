

%-----------------
% HOUSE KEEPING
%----------------

clear all;  
clc; 
close all;

addpath('./auxfiles')
addpath('./results')
acpt =0
%------------------------------------------------------------
% SETTINGS
%------------------------------------------------------------
rng(2666);;

model_vec = {'m2lev'};     % Select Model
instrList = {'EGON_KUTTNER_NI'};        % {'EGON_KUTTNER'}
p = 12;                                 % Number of lags
nex_ = 1;                               % Constant
T0 = p+36;                              % Length of pre-sample for Minnesota Prior
str_sample_init = '1990-01-01';         % Starting date of the sample (include pre-sample)
str_sample_end  = '2007-06-01';         % End date of the sample
str_iv_init     = '1994-01-01';         % Starting date of the sample for the proxy
MP = 0;                                 % Tightness of the Minnesota Prior (0: loose; 1: tight)
i_var_instr = instrList;
i_var_str =  {'FFR_SSR', 'IPM','UNRATE', 'PPI_FIN', 'BAA_10YMOODY'};
i_var_transf =  {};     % Variable that require additional transformations
nCalc = length(i_var_transf);
i_var_str_names =  i_var_str; % Name of variables (for plots)
varSelec = [1 2 3 4 5]; % Select variables to plot

nlags_ = p;

data_file = 'vardata_extended';
data_spreadsheet = 'Sheet1';
mmodel = model_vec(:,1);
%-------------------------------------------
% Load data and construst dummy observations
%-------------------------------------------
vm_dummy;

Mobs       = size(mm,1);            
MM         = [ones(Mobs,1) mm]; % Add constant to matrix of proxies

% For uninformative prior uncomment the following two lines
% XXdum = [];
% YYdum = [];

%-------------------------------------------
% Declare objects for estimation
%-------------------------------------------

n = nv;
nIV = size(i_var_instr,2);
nshocks=nIV;

e = eye(n); % create identity matrix
aalpha_index = 2:n;
ddelta_index = 1;
a = cell(p,1);

% Define matrices to compute IRFs      
J = [eye(n);repmat(zeros(n),p-1,1)]; % Page 12 RWZ
F = zeros(n*p,n*p);    % Matrix for Companion Form
I  = eye(n);
for i=1:p-1
    F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
end

% Estimation Preliminaries
X = [XXdum; XXact];
Y = [YYdum; YYact];
T = size(X, 1);
ndum = size(XXdum, 1);
nex = nex_;

% Compute OLS estimates
B = (X'*X)\(X'*Y); % Point estimates
U = Y-X*B;         % Residuals
Sigmau = U'*U/(T-p*n-1);   % Covariance matrix of residuals
F(1:n,1:n*p)    = B(1:n*p,:)';

%------------------------------------------------------------
% MCMC Algorithm
%------------------------------------------------------------

% set preliminaries for priors
N0=zeros(size(X',1),size(X,2));
nnu0=0;
nnuT = T +nnu0;
NT = N0 + X'*X;    
Bbar0=B;
S0=Sigmau;
BbarT = NT\(N0*Bbar0 + (X'*X)*B);
ST = (nnu0/nnuT)*S0 + (T/nnuT)*Sigmau + (1/nnuT)*((B-Bbar0)')*N0*(NT\eye(n*p+nex))*(X'*X)*(B-Bbar0); %% Constant (check)
STinv = ST\eye(n);

record=0;     
counter = 0;
fflagEXP = 0;

% Drop constant from M
MM = MM(:, 2:end);

bet = 0.01;
signu = 0.04;

[Q, ~] = qr(randn(n, nIV));
Xstar = randn(n, 1);
%Q = Xstar / norm(Xstar);

R=mvnrnd(zeros(n,1),STinv/nnuT,nnuT)';
Sigmadraw=(R*R')\eye(n);
bbeta = B(:);
SigmaB = kron(Sigmadraw,NT\eye(n*p+nex_));
SigmaB = (SigmaB+SigmaB')/2;
Bdraw = mvnrnd(bbeta,SigmaB);

Bdraw= reshape(Bdraw,n*p+nex,n); % Reshape Bdraw from vector to matrix
Udraw = Y-X*Bdraw;      % Store residuals for IV regressions
LC =chol(Sigmadraw,'lower');
A0chol = (LC')\eye(size(LC,1));

lnp0 = loglik_m_given_y(MM, Udraw(ndum+1:end, :), LC, Q(:, 1), bet, signu)

A0 = (A0chol*Q);
Ap = Bdraw*A0;;
A0
paravec = [A0(:); Ap(:); bet];

save('-ascii', 'test_para.txt', 'paravec')
save('-ascii', 'test_lik.txt', 'lnp0')

