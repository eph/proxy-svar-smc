tic

%-----------------
% HOUSE KEEPING
%----------------

clear all;  
clc; 
close all;

addpath('./auxfiles')
addpath('./results')

%------------------------------------------------------------
% SETTINGS
%------------------------------------------------------------
model_vec = {'4eq'};                    % Select Model ('4eq' or '5eq')
instrList = {'MHF'};
p = 12;                                 % Number of lags
nex_ = 1;                               % Constant
Horizon = 48;                           % Horizon for calculation of impulse responses
MP = 0;                                 % Tightness of the Minnesota Prior (0: loose; 1: tight)
T0 = p+36;                              % Length of pre-sample for Minnesota Prior
nd = 100000;                            % Number of draws in MC chain
bburn = 0.1*nd;                         % Burn-in period
fflagFEVD = 0;                          % Compute FEVD
printFig = 1;                           % Save Figures in PDF

str_sample_init = '1990-01-01';         % Starting date of the sample (include pre-sample)
str_sample_end  = '2007-06-01';         % End date of the sample
str_iv_init     = '1994-01-01';         % Starting date of the sample for the proxy


ptileVEC = [0.05 0.16 0.50 0.86 0.95]; % Percentiles

%------------------------------------------------------------
% PRIOR SELECTION -- SIGMA_NU
% IG - Inverse Gamma
% FI - Truncated at pr_trunc x std(M_t)
%------------------------------------------------------------
prior_type = 'IG';
pr_trunc = 0.5;                      % only used for FIXED

mu0 = 0.00;                             % prior mean
V0 = 0.1^2;                             % prior variance

s0 = 0.02;
nu0 = 2.00;
%------------------------------------------------------------
% Sampler Settings
%------------------------------------------------------------
% this is the mixture probability for the "RW" IG propsal on SIGMA
rwmh_sigma_prob = 0.0;
rwmh_df = 5;                            % Tune-up parameter for mixture proposal distribution for (\Phi,\Sigma)
nu = rwmh_df;                           % Tune-up parameter for mixture proposal distribution for (\Phi,\Sigma)

%----------------------------------------------------------------------
% LOAD DATA
%---------------------------------------------------------------------
for mCounter = 1:size(model_vec,2)
    mmodel = model_vec(:,mCounter);
    if strcmp(mmodel, '5eq')
        i_var_instr = instrList;
        i_var_str =  {'EFFR_LW', 'LIPM', 'UNRATE', 'LPPI', 'BAA10YMOODY'};
        i_var_transf =  {};     % Variable that require additional transformations
        nCalc = length(i_var_transf);
        i_var_str_names =  i_var_str; % Name of variables (for plots)
        varSelec = [1 2 3 4 5]; % Select variables to plot
    elseif strcmp(mmodel, '4eq')
        i_var_instr = instrList;
        i_var_str =  {'EFFR_LW', 'LIPM', 'UNRATE', 'LPPI'};
        i_var_transf =  {};     % Variable that require additional transformations
        nCalc = length(i_var_transf);
        i_var_str_names =  i_var_str; % Name of variables (for plots)
        varSelec = [1 2 3 4]; % Select variables to plot
    end

    nlags_ = p;

    data_file = '../data/CHdata.txt';
    data_spreadsheet = 'Sheet1';

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
    nnuT = T - n*p -nex_;%# +nnu0;
    NT = N0 + X'*X;    
    Bbar0=B;
    S0=Sigmau;
    BbarT = NT\(N0*Bbar0 + (X'*X)*B);
    ST = (nnu0/nnuT)*S0 + (T/nnuT)*Sigmau + (1/nnuT)*((B-Bbar0)')*N0*(NT\eye(n*p+nex))*(X'*X)*(B-Bbar0); %% Constant (check)
    STinv = ST\eye(n);

    record=0;     
    counter = 0;
    fflagEXP = 0;

    disp('                                                                  ');
    disp('        BAYESIAN ESTIMATION OF VAR VIA BLOCK MCMC                 ');
    disp('                                                                  ');

    % Drop constant from M
    MM = MM(:, 2:end);

    bet = 0.01;
    signu = 0.04;

    %     [Q, ~] = qr(randn(n, nIV));
    Xstar = randn(n, 1);
    Q = Xstar / norm(Xstar);
    
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



    lnp0 = loglik_m_given_y(MM, Udraw(ndum+1:end, :), LC, Q, bet, signu);

    % Initialize Omega1 for IRFs
    Omega1 = [LC(:,nshocks);zeros((p-1)*n,nshocks)];

    nq = 1;
    acpt_rf   = 0;
    acpt_Q    = 0;


    Fstar = F;
    Sigmadrawstar = Sigmadraw;


    % MCMC Chain 
    % Define objects that store the draws
    Ltilde = zeros(nd-bburn,Horizon+1,n,nshocks);                      % define array to store IRF
    LtildeAdd = zeros(nd-bburn,Horizon+1,n+nCalc,nshocks);
    irfCalc = zeros(nd-bburn,Horizon+1,nCalc,nshocks);                     % store labor productivity IRF
    W = zeros(nd-bburn,Horizon+1,n+nCalc,nshocks);                         % define array to store FVD
    EETA  = zeros(nd-bburn,n-1);
    ppsi_levels = zeros(nd-bburn,n);
    ppsi_diff   = zeros(nd-bburn,n-1);
    REL = zeros(nd-bburn,1);
    BET = zeros(nd-bburn,1);
    SIG = zeros(nd-bburn,1);
    A0MAT = zeros(nd-bburn, n, n);
    ApMAT = zeros(nd-bburn, n*p+nex, n);
    while record<nd

        %------------------------------------------------------------
        % Gibbs Sampling Algorithm
        %------------------------------------------------------------
        % STEP ONE: Draw from the B, SigmaB | Y
        %------------------------------------------------------------
        % Step 1: Draw from the marginal posterior for Sigmau p(Sigmau|Y,X)
        if (rand()<rwmh_sigma_prob)     % draw sigma* ~ IW(sigma, nu), phi*|sigma*, Y
            nu = nobs;
            Sigmadrawstar = iwishrnd(nu*Sigmadraw, nu+n+1);

            bbeta = B(:);
            SigmaB = kron(Sigmadrawstar,NT\eye(n*p+nex));
            SigmaB = (SigmaB+SigmaB')/2;
            Bdrawstar = mvnrnd(bbeta,SigmaB);
            
            Bdrawstar= reshape(Bdrawstar,n*p+nex,n);% Reshape Bdraw from vector to matrix
            Udrawstar = Y-X*Bdrawstar;              % Store residuals for IV regressions
            
            LCstar    = chol(Sigmadrawstar,'lower');
            A0cholstar    = (LCstar')\eye(size(LCstar,1));
            
            Fstar(1:n,1:n*p)    = Bdrawstar(1:n*p,:)';

            lnp1 = loglik_m_given_y(MM, Udrawstar(ndum+1:end, :), LCstar, Q(:,nIV), bet, signu);

            q1 = pdf_ln_iwish(nu*Sigmadraw, nu+n+1, Sigmadrawstar);
            q0 = pdf_ln_iwish(nu*Sigmadrawstar, nu+n+1, Sigmadraw);

            psigma1 = pdf_ln_iwish(ST, nnuT, Sigmadrawstar);
            psigma0 = pdf_ln_iwish(ST, nnuT, Sigmadraw);            

            % metropolis hastings
            alp = exp((lnp1+psigma1) - (lnp0+psigma0) - (q1 - q0));


        else                            % draw phi*, sigma* ~ phi, sigma | Y

            R=mvnrnd(zeros(n,1),STinv/nnuT,nnuT)';
            Sigmadrawstar=(R*R')\eye(n);

            % Step 2: Taking newSigma as given draw for B using a multivariate normal    
            bbeta = B(:);
            SigmaB = kron(Sigmadrawstar,NT\eye(n*p+nex));
            SigmaB = (SigmaB+SigmaB')/2;
            Bdrawstar = mvnrnd(bbeta,SigmaB);
            % Storing unrestricted draws
            
            Bdrawstar= reshape(Bdrawstar,n*p+nex,n);% Reshape Bdraw from vector to matrix
            Udrawstar = Y-X*Bdrawstar;      % Store residuals for IV regressions
            
            LCstar    = chol(Sigmadrawstar,'lower');
            A0cholstar    = (LCstar')\eye(size(LCstar,1));
            
            Fstar(1:n,1:n*p)    = Bdrawstar(1:n*p,:)';

            lnp1 = loglik_m_given_y(MM, Udrawstar(ndum+1:end, :), LCstar, Q(:,nIV), bet, signu);

            % metropolis hastings
            alp = exp(lnp1 - lnp0);

        end 

        if rand < alp
            lnp0 = lnp1;
            
            Udraw = Udrawstar;
            LC = LCstar;
            A0chol = A0cholstar;
            
            F = Fstar;
            Sigmadraw = Sigmadrawstar;
            Bdraw = Bdrawstar;           
            acpt_rf = acpt_rf+1;
        end
        
        %------------------------------------------------------------
        % STEP TWO: Draw from Q distribution
        %------------------------------------------------------------
        %            [Qstar, ~] = qr(randn(n, nIV));
        Xstar = randn(n, 1);
        Qstar = Xstar / norm(Xstar);
        
        lnp1 = loglik_m_given_y(MM, Udraw(ndum+1:end, :), LC, Qstar(:,nIV), bet, signu);

        % metropolis hastings
        alp = exp(lnp1 - lnp0);
        if rand < alp
            lnp0 = lnp1;
            Q = Qstar;
            acpt_Q = acpt_Q+1;
        end

        % normalize A0
        A0 = A0chol*Q;
        if A0(1) < 0
            Q = -Q;
        end

        %------------------------------------------------------------
        % STEP THREE: BET and SIGMA from MVN-IG / MVN-IW
        %------------------------------------------------------------
        scale_mat = Q' / LC;
        Xe = (scale_mat*Udraw(ndum+1:end, :)')';
 
        % bayesian linear regresion model
        Bhat = inv(Xe'*Xe)*Xe'*MM;
        Vp = inv(Xe'*Xe + inv(V0/signu^2));
        mup = Vp * (Xe'*Xe*Bhat + inv(V0)*mu0);

        if (prior_type == 'IG')          
            nu1 = size(Xe, 1)-1 + nu0;
            s1 = ( (MM-Xe*mup)'*(MM-Xe*mup) + nu0*s0^2)/nu1;

            signu = igrand(s1, nu1);
        elseif (prior_type == 'FI')     
            signu = pr_trunc*std(MM);
        end
 
        bet = mvnrnd(mup, signu^2*Vp);


        lnp0 = loglik_m_given_y(MM, Udrawstar(ndum+1:end, :), LCstar, Q(:,nIV), bet, signu);

        record=record+1;
        counter = counter +1;
        if counter==0.05*nd
            disp(['         DRAW NUMBER:   ', num2str(record)]);
            disp('                                                                  ');
            disp(['     REMAINING DRAWS:   ', num2str(nd-record)]);
            disp('                                                                  ');
            counter = 0;
            
            fprintf('RF acpt rate: %5.3f\n',   acpt_rf/record)
            fprintf('Q acpt rate: %5.3f\n',    acpt_Q/(record*nq))
        end

        if record > bburn
            
            ffactor = LC*Q;
            A0 = A0chol*Q;
            Aplus = Bdraw(1:n*p,:)*A0;


            % Compute cumulative coefficients 
            a0 = A0(:,1);
            for l=1:p
                a{l} = Aplus((l-1)*n+1:l*n,1);
            end
            aalpha = zeros(p+1,n-1);
            for j=1:n-1
                jj=aalpha_index(1,j);
                aalpha(1,j)  = -((e(:,jj)'*a0));
                for l=1:p
                    aalpha(l+1,j)  = ((e(:,jj)'*a{l}));
                end
            end
            ddelta = zeros(p+1,1);
            ddelta(1,1) = ((e(:,ddelta_index )'*a0));
            for l=1:p
                ddelta(l+1,1)  = -((e(:,ddelta_index )'*a{l}));
            end
            tmp_levels   = zeros(1,n-1);

            for l=0:p
                tmp_levels = aalpha(l+1,1:n-1) + tmp_levels;
            end
            % variables in levels
            tmp_diff   = zeros(1,n-1);
            for l=0:p
                for ii=0:l
                    tmp_diff = aalpha(ii+1,:) +  tmp_diff ;
                end
            end
            tmpR=0;
            for l=0:p
                tmpR=ddelta(l+1,1)+tmpR;
            end

            aalphaEta = aalpha./a0(1);
            ddeltaEta = ddelta./a0(1);
            den = ddeltaEta(1) + sum(ddeltaEta(2:end));       
            num = sum(aalphaEta);
            numtemp = sum(cumsum(aalphaEta));
        
            REL(record-bburn,1) = bet^2/(bet^2 + signu^2);
            EETA(record-bburn,:) = -A0(2:n,1)./A0(1,1);
            ppsi_levels(record-bburn,:) = [sum(ddeltaEta(2:end)) num];
            ppsi_diff(record-bburn,1:n-1)   = numtemp;
            BET(record-bburn,1) = bet;
            SIG(record-bburn,1) = signu;
            IRF_T    = vm_irf(F,J,ffactor,Horizon+1,n,Omega1);        
            Ltilde(record-bburn,:,:,:) = IRF_T(1:Horizon+1,:,:);

            if nCalc
                for ii = 1:length(i_transf)
                    irfCalc(record-bburn,:,ii,:) = cumsum(squeeze(IRF_T(:,i_transf(ii),:)));
                end
            end

            if fflagFEVD ==1
                W(record-bburn,:,:,:)=variancedecompositionFD(F,J,Sigmadraw,ffactor,n,Horizon,i_transf);
            end
        end
    end
    
    
    LtildeAdd(:,:,1:n,:) = Ltilde;
    LtildeAdd(:,:,n+1:n+nCalc,:) = irfCalc;
    
    SVAR.LtildeImpact = LtildeAdd(:,1,:,:);
   
    LtildeFull = quantile(LtildeAdd,ptileVEC);
    SVAR.LtildeFull = permute(LtildeFull,[3,2,1,4]);

    WhFull = quantile(W,ptileVEC);
    SVAR.WhFull = permute(WhFull,[3,2,1,4]);
      
    SVAR.EETAFull = quantile(EETA,ptileVEC);
    SVAR.PPSIFull = quantile(ppsi_levels,ptileVEC);
    SVAR.PPSIDFull = quantile(ppsi_diff,ptileVEC);
    SVAR.RELFull  = quantile(REL,ptileVEC);
    SVAR.BETFull  = quantile(BET,ptileVEC);
    SVAR.SIGFull  = quantile(SIG,ptileVEC);
    SVAR.SIG = SIG;
    SVAR.BET = BET;
    SVAR.acpt_rf = acpt_rf/record;
    SVAR.acpt_Q  = acpt_Q/(record*nq);
    SVAR.p = p;
    SVAR.nd = nd;
    SVAR.bburn = bburn;
    SVAR.i_var_str_names = i_var_str_names;
    SVAR.fflagFEVD = fflagFEVD;
    SVAR.mmodel = mmodel;
    SVAR.pr_trunc = 0
    SVAR.pr_truncFlag = 0;


    SVAR.printFig = printFig;
    SVAR.i_var_instr = i_var_instr;
    SVAR.EETA = EETA;
    SVAR.varSelec = varSelec;
    
    savefileirf = strcat('./results/Result_',char(mmodel),char(i_var_instr),'p_',num2str(p),'_pr_',prior_type,'MP_',num2str(MP),'.mat');
    save(savefileirf,'SVAR'); 
    toc

    vm_plot_irf_bvar(mmodel,prior_type,i_var_instr,p,MP,fflagFEVD,printFig)
    diary off
end
