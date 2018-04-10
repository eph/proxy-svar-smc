tic


%-----------------
% Housekeeping    
%----------------

clear all;  
close all;
clc; 

addpath('results')
addpath('auxfiles')
printFig = 0;
model_vec = {'Chol_4eq'};
nd = 5000;        % Number of draws in MC chain         
instrList = {'MHF'};
prior_type = '';
    
for iCounter = 1:size(instrList,2)
    for mCounter = 1:size(model_vec,2)
        mmodel = model_vec(:,mCounter);
       if strcmp(mmodel, 'Chol_4eq')
            i_var_instr = instrList(:,iCounter);
            i_var_str =  {'EFFR_LW', 'LIPM', 'UNRATE', 'EFFR'};
            i_var_transf =  {};
            nCalc = length(i_var_transf);
            Horizon = 48;
            i_var_str_names =  i_var_str;
            str_sample_init = '1990-01-01';
            str_sample_end  = '2007-06-01';
            str_iv_init     = '1994-01-01';
            varSelec = [1 2 3 4];
            nIV = size(i_var_instr,2);
            nshocks=length(i_var_str); % Number of shocks to identify
            p = 12;                     % Number of lags
            aalpha_index = 1:3;
            ddelta_index = 4; % ffr
       elseif strcmp(mmodel, 'Chol_5eq')
            i_var_instr = instrList(:,iCounter);
            i_var_str =  {'EFFR_LW', 'LIPM', 'UNRATE', 'BAA10YMOODY', 'EFFR'}; %,'BAA_10YMOODY'
            i_var_transf =  {};
            nCalc = length(i_var_transf);
            Horizon = 48;
            i_var_str_names =  i_var_str;
            str_sample_init = '1990-01-01';
            str_sample_end  = '2007-06-01';
            str_iv_init     = '1994-01-01';
            varSelec = [1 2 3 4 5];
            nIV = size(i_var_instr,2);
            nshocks=length(i_var_str); % Number of shocks to identify
            p = 12;                     % Number of lags
            aalpha_index = 1:4;
            ddelta_index = 5; % ffr
       elseif strcmp(mmodel, 'Coibion')
            i_var_instr = instrList(:,iCounter);
            i_var_str =  {'LIPM', 'UNRATE', 'BAA10YMOODY', 'EFFR_LW','MRR'};
            i_var_transf =  {};
            nCalc = length(i_var_transf);
            Horizon = 48;
            i_var_str_names =  i_var_str
            str_sample_init = '1994-01-01';
            str_sample_end  = '2007-06-01';
            str_iv_init     = '1994-01-01';
            varSelec = [1 2 3 4 5];
            nIV = size(i_var_instr,2);
            nshocks=length(i_var_str); % Number of shocks to identify
            p = 12;                     % Number of lags
            aalpha_index = 1:4;
            ddelta_index = 5; % ffr
       elseif strcmp(mmodel, 'CoibionNewInstr')
            i_var_instr = instrList(:,iCounter);
            i_var_str =  {'IPM','PPI_FIN','UNRATE','CRB_COMMP','RR_SHOCK_SUB'};
            i_var_transf =  {};
            nCalc = length(i_var_transf);
            Horizon = 48;
            i_var_str_names =  {'IPM','PPI_FIN','UNRATE','CRB_COMMP','RR_SHOCK_SUB'};
            str_sample_init = '1994-01-01';
            str_sample_end  = '2007-06-01';
            str_iv_init     = '1994-01-01';
            varSelec = [1 2 3 4 5];
            nIV = size(i_var_instr,2);
            nshocks=length(i_var_str); % Number of shocks to identify
            p = 12;                     % Number of lags
            aalpha_index = 1:4;
            ddelta_index = 5; % ffr
       elseif strcmp(mmodel, 'CoibionWithSpread')
            i_var_instr = instrList(:,iCounter);
            i_var_str =  {'IPM','PPI_FIN','UNRATE','CRB_COMMP','BAA_10YMOODY','RR_SHOCK_NOSPREAD'};
            i_var_transf =  {};
            nCalc = length(i_var_transf);
            Horizon = 48;
            i_var_str_names =  {'IPM','PPI_FIN','UNRATE','CRB_COMMP','BAA_10YMOODY','RR_SHOCK_NOSPREAD'};
            str_sample_init = '1994-01-01';
            str_sample_end  = '2007-06-01';
            str_iv_init     = '1994-01-01';
            varSelec = [1 2 3 4 5 6];
            nIV = size(i_var_instr,2);
            nshocks=length(i_var_str); % Number of shocks to identify
            p = 12;                     % Number of lags
            aalpha_index = 1:5;
            ddelta_index = 6; % ffr
        elseif strcmp(mmodel, 'CoibionNewInstrWithSpread')
            i_var_instr = instrList(:,iCounter);
%             i_var_str =  {'IPM','PPI_FIN','UNRATE','CRB_COMMP','RR_SHOCK_SUB','BAA_10YMOODY'};
%             i_var_str_names =  {'IPM','PPI_FIN','UNRATE','CRB_COMMP','RR_SHOCK_SUB','BAA_10YMOODY'};
            i_var_str =  {'IPM','PPI_FIN','UNRATE','CRB_COMMP','BAA_10YMOODY','RR_SHOCK_SUB'};
            i_var_str_names =  {'IPM','PPI_FIN','UNRATE','CRB_COMMP','BAA_10YMOODY','RR_SHOCK_SUB'};
            i_var_transf =  {};
            nCalc = length(i_var_transf);
            Horizon = 48;
            
            str_sample_init = '1994-01-01';
            str_sample_end  = '2007-06-01';
            str_iv_init     = '1994-01-01';
            varSelec = [1 2 3 4 5 6];
            nIV = size(i_var_instr,2);
            nshocks=length(i_var_str); % Number of shocks to identify
            p = 12;                     % Number of lags
%             aalpha_index = [1 2 3 4 6];
%             ddelta_index = 5; % ffr
            aalpha_index = [1 2 3 4 5];
            ddelta_index = 6; % ffr
        end

    data_file = '../data/CHdata.txt';
    data_spreadsheet = 'Sheet1';

    nlags_ = p;
    T0 = p+36;
    nex_ = 1;
    MP=0;
    
    vm_dummy;
    
    n = size(i_var_str,2);  % Number of Endogenous variables
    nobs = size(YY,1)-T0;
    e=eye(n);
    bburn = 0.2*nd;
    ptileVEC = [0.05 0.16 0.50 0.86 0.95];
    fflagFEVD = 0;


    Mobs       = size(mm,1);            
    MM         = [ones(Mobs,1) mm];
    
    

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
    U = Y-X*B;      % Residuals
    Sigmau = U'*U/(T-p*n-1);   % Covariance matrix of residuals
    LC =chol(Sigmau,'lower');
    A0chol = (LC')\eye(size(LC,1));
    F(1:n,1:n*p)    = B(1:n*p,:)';

    %------------------------------------------------------------
    % MCMC settings
    %------------------------------------------------------------
    bburn = 0.2*nd;
    ptileVEC = [0.05 0.16 0.50 0.86 0.95];

    % set preliminaries for priors
    N0=zeros(size(X',1),size(X,2));
    nnu0=0;
    nnuT = T - n*p -nex_;%# +nnu0; %T +nnu0;
    NT = N0 + X'*X;    
    Bbar0=B;
    S0=Sigmau;
    BbarT = NT\(N0*Bbar0 + (X'*X)*B);
    ST = (nnu0/nnuT)*S0 + (T/nnuT)*Sigmau + (1/nnuT)*((B-Bbar0)')*N0*(NT\eye(n*p+nex))*(X'*X)*(B-Bbar0); %% Constant (check)
    STinv = ST\eye(n);


    % MCMC Chain 

    % initialize Omega1 for IRFs
    Omega1 = [LC;zeros((p-1)*n,size(LC,2))];
    % Define objects that store the draws

    Ltilde = zeros(nd-bburn,Horizon+1,n,nshocks);                      % define array to store IRF
    LtildeAdd = zeros(nd-bburn,Horizon+1,n+nCalc,nshocks);
    irfCalc = zeros(nd-bburn,Horizon+1,nCalc,nshocks);                     % store labor productivity IRF
    W = zeros(nd-bburn,Horizon+1,n+nCalc,nshocks);                         % define array to store FVD
    EETA = zeros(nd-bburn,n);
    ppsi_levels = zeros(nd-bburn,n);
    ppsi_diff   = zeros(nd-bburn,n-1);

    
    
    a = cell(p,1);
    record=0;     
    counter = 0;
    fflagEXP = 0;
    
    disp('                                                                  ');
    disp('        BAYESIAN ESTIMATION OF VAR: DIRECT SAMPLING...            ');
    disp('                                                                  ');

    while record<nd

        % Gibbs Sampling Algorithm

        % Step 1: Draw from the marginal posterior for Sigmau p(Sigmau|Y,X)
        R=mvnrnd(zeros(n,1),STinv/nnuT,nnuT)';
        Sigmadraw=(R*R')\eye(n);

        % Step 2: Taking newSigma as given draw for B using a multivariate normal    
        bbeta = B(:);
        SigmaB = kron(Sigmadraw,NT\eye(n*p+nex));
        SigmaB = (SigmaB+SigmaB')/2;
        Bdraw = mvnrnd(bbeta,SigmaB);

        % Storing unrestricted draws
        Bdraw= reshape(Bdraw,n*p+nex,n); % Reshape Bdraw from vector to matrix
        Udraw = Y-X*Bdraw;      % Store residuals for IV regressions
        LC =chol(Sigmadraw,'lower');
        A0    = (LC')\eye(size(LC,1));
        Aplus = Bdraw(1:n*p,:)*A0;
        F(1:n,1:n*p)    = Bdraw(1:n*p,:)';
        
        eigen           = eig(F);
        eigen           = max(eigen);
        largeeig       = abs(eigen);
        if largeeig > 1
%             continue
            fflagEXP = fflagEXP + 1;
        end

        ffactor  = LC;
        
        % Reliability (not sure how to compute PhiPhip)
        Sigmm   = mm'*mm/Mobs;
        ED      = eye(nIV)*sum(sum(mm,2)~=0)/Mobs;
        
        a0 = A0(:,ddelta_index);
        for l=1:p
            a{l} = Aplus((l-1)*n+1:l*n,ddelta_index);
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
        % in levels
        for l=0:p
            tmp_levels = aalpha(l+1,1:n-1) + tmp_levels;
        end
    % in first differences
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
        
        aalphaEta = aalpha./a0(ddelta_index);
        ddeltaEta = ddelta./a0(ddelta_index);
        den = ddeltaEta(1) + sum(ddeltaEta(2:end));       
        num = sum(aalphaEta);
        numtemp = sum(cumsum(aalphaEta));

        record=record+1;
        counter = counter +1;
        if counter==0.05*nd
            disp(['         DRAW NUMBER:   ', num2str(record)]);
            disp('                                                                  ');
            disp(['     REMAINING DRAWS:   ', num2str(nd-record)]);
            disp('                                                                  ');
            counter = 0;
        end
        if record > bburn
            if strcmp(mmodel, 'Chol_5eq')
            A0temp = A0;    
            A0temp(4,:) = A0(4,:)/A0(4,4);
            ffactor = inv(A0temp');
            ffactor(:,4) = ffactor(:,4)*0.1;
            end
            
            IRF_T    = vm_irf(F,J,ffactor,Horizon+1,n,Omega1);
            Ltilde(record-bburn,:,:,:) = IRF_T(1:Horizon+1,:,:);
            EETA(record-bburn,:) = -A0(:,ddelta_index)./A0(ddelta_index,ddelta_index);
            ppsi_levels(record-bburn,:) = [sum(ddeltaEta(2:end)) num];
            ppsi_diff(record-bburn,1:n-1)   = numtemp;
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

    SVAR.EETAFull = quantile(EETA,ptileVEC);
    SVAR.PPSIFull = quantile(ppsi_levels,ptileVEC);
    SVAR.PPSIDFull = quantile(ppsi_diff,ptileVEC);
     
    WhFull = quantile(W,ptileVEC);
    SVAR.WhFull = permute(WhFull,[3,2,1,4]);
    
    SVAR.varSelec = varSelec;
    SVAR.i_var_str_names = i_var_str_names;
    SVAR.EETA = EETA;
    
    savefileirf = strcat('./results/Result_',char(mmodel),char(i_var_instr),'p_',num2str(p),'_pr_',prior_type,'MP_',num2str(MP),'.mat');
    save(savefileirf,'SVAR'); 
    toc

    vm_plot_irf_bvar(mmodel,prior_type,i_var_instr,p,MP,fflagFEVD,printFig)

    end
end
