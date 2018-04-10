%==========================================================================
% This script produces Figure 7 in Caldara, Dario and Edward Herbst (2018),
% "Monetary Policy, Real Activity, and Credit Spreads:
% Evidence from Bayesian Proxy SVAR", American Economic Journal: Macroeconomics.
% The code is take from the paper "Proxy SVARs: Asymptotic Theory, 
% Bootstrap Inference, and the Effects of Income Tax Changes in the United 
% States" by Carsten Jentsch and Kurt G. Lunsford.  
% 
% This script calls the following workspaces:
%   VARData.mat 
%   ProxyData.mat
% 
% This script call the following functions:
%   IdentifyFunc.m
%   acf.m
%==========================================================================
close all
clc
clear;
format short g
tic;

%==========================================================================
% Parameters of the Program
%==========================================================================
%model parameters
lags = 12;              %number of lags in the structural VAR
nImp = 49;              %number of steps in the IRF
alpha = 0.9;            %level of significance for confidence intervals

%bootstrap parameters
nBoot = 20000;          %number of bootstrap replications
BlockSize = 10;         %size of blocks in the bootstrap
seed = 1;               %seed for random number generator
rng(seed);              %iniate the random number generator

%program toggles
toggleBoot = 2;         %1 to use moving block bootstrap
                        %2 to use wild bootstrap (Mertens Ravn AER)

varSize = 4;            %select number of variables (4 or 5)
%==========================================================================
% Load Data
%==========================================================================
data_file = '../data/CHdata.txt';
data_spreadsheet = 'Sheet1';
i_var_transf =  {};

str_sample_init = '1994-01-01';
str_sample_end  = '2007-06-01';
str_iv_init     = '1994-01-01';

i_var_instr = {'MHF'};

if varSize == 4
    i_var_str =  {'EFFR_LW', 'LIPM','UNRATE','LPPI'};
elseif varSize == 5
    i_var_str =  {'EFFR_LW', 'LIPM','UNRATE','LPPI','BAA10YMOODY'};
end

vm_loaddata

%dimensions of the VAR data
[T,K] = size(YY);

M = mm(1+lags:T,:);  %adjust the proxy data for the number of lags

T0 = T-size(M,1)-lags;
r = size(M,2);              %number of proxy variables

%==========================================================================
% Estimate the VAR
%==========================================================================

%left-hand side variables
LHS = YY(1+lags:T,:);

%right-hand side variables
preRHS = ones(T-lags,1);
for j = 1:lags
    preRHS = [preRHS,YY(1+lags-j:T-j,:)];
end
RHS = preRHS;

%VAR coefficients and innovations
A_hat = (RHS'*RHS)\(RHS'*LHS);
U_hat = LHS - RHS*A_hat;


%==========================================================================
% Autocorrelation Functions of VAR Residuals
%==========================================================================
%acf of u(t)
rhoU = acf(U_hat,7);

%acf of |u(t)|
rhoUabs = acf(abs(U_hat),7);

%acf of u(t)^2
rhoUsq = acf((U_hat.^2),7);

%table of autocorrelations for appendix
Table = [rhoU;rhoUabs;rhoUsq];


%==========================================================================
% Identify the Structural Shocks
%==========================================================================
%construct B1 for both the APITR and ACITR ordered first
U_hat = U_hat(T0+1:end,:,:);
H1_hat = zeros(K,r,2);

H1_hat = IdentifyFunc(U_hat,M);


%==========================================================================
% Produce the IRFs
%==========================================================================
IRF = zeros(K,nImp);

%generate the IRFs
IRF(1:K,1) = H1_hat;
history = [IRF(1:K,1);zeros((lags-1)*K,1)];
for j = 2:nImp
    IRF(1:K,j) = A_hat(2:K*lags+1,:)'*history;
    history = [IRF(1:K,j);history(1:(lags-1)*K,1)];
end



%==========================================================================
% Bootstrap and Confidence Intervals
%==========================================================================
%set the number of blocks
if toggleBoot == 1
    nBlock = ceil((T-lags)/BlockSize);
end

%arrays to store relevant variables
H1_hatboot = zeros(K,r,nBoot,2);
IRF_boot = zeros(K,nImp,nBoot,2);
M_count = zeros(nBoot,r,2);

for order = 1 %original code allowed for multiple instruments and orderings
    
    %if using the moving block bootstrap, create the blocks and centerings
    if toggleBoot == 1
        Blocks = zeros(BlockSize,K,T-lags-BlockSize+1);
        MBlocks = zeros(BlockSize,r,T-lags-BlockSize+1);
        for j = 1:T-lags-BlockSize+1-T0
            Blocks(:,:,j) = U_hat(j:BlockSize+j-1,:,order);
            MBlocks(:,:,j) = M(j:BlockSize+j-1,:);
        end

        %center the bootstrapped VAR errors
        centering = zeros(BlockSize,K);
        for j = 1:BlockSize
            centering(j,:) = mean(U_hat(j:T-lags-BlockSize+j-T0,:,order),1);
        end
        centering = repmat(centering,[nBlock,1]);
        centering = centering(1:T-lags-T0,:);

        %center the bootstrapped proxy variables
        Mcentering = zeros(BlockSize,r);
        for j = 1:BlockSize
            subM = M(j:T-lags-BlockSize+j-T0,:);
            Mcentering(j,:) = mean(subM((subM(:,1) ~= 0),1),1);
        end
        Mcentering = repmat(Mcentering,[nBlock,1]);
        Mcentering = Mcentering(1:T-lags-T0,:);
    end

    %loop for the bootstrap replications
    for boot = 1:nBoot
        if toggleBoot == 1
            %draw bootstrapped residuals and proxies
            index = ceil((T - lags - BlockSize + 1)*rand(nBlock,1));
            U_boot = zeros(nBlock*BlockSize,K);
            M_boot = zeros(nBlock*BlockSize,r);
            for j = 1:nBlock
                U_boot(1+BlockSize*(j-1):BlockSize*j,:) = Blocks(:,:,index(j,1));
                M_boot(1+BlockSize*(j-1):BlockSize*j,:) = MBlocks(:,:,index(j,1));
            end
            U_boot = U_boot(1:T-lags-T0,:);
            M_boot = M_boot(1:T-lags-T0,:);

            %center the bootstrapped residuals and proxies
            U_boot = U_boot - centering;
            for j = 1:r
                M_boot((M_boot(:,j)~=0),j) =...
                    M_boot((M_boot(:,j)~=0),j) - Mcentering((M_boot(:,j)~=0),j);
            end

            %produce the bootstrapped left- and right-hand side data
            LHS_boot = zeros(T-lags-T0,K);
            RHS_boot = [ones(T-lags-T0,1),zeros(T-lags-T0,K*lags)];
            RHS_boot(1,:) = RHS(1,:,order);
            for j = 1:T-lags-1-T0
                LHS_boot(j,:) = RHS_boot(j,:)*A_hat(:,:,order) + U_boot(j,:);
                RHS_boot(j+1,:) = [1,LHS_boot(j,:),RHS_boot(j,2:K*(lags-1)+1)];
            end
            LHS_boot(T-lags-T0,:) =...
                RHS_boot(T-lags-T0,:)*A_hat(:,:,order) + U_boot(T-lags-T0,:);
        elseif toggleBoot == 2
            %generate the shocks for the wild bootstrap
            E = 2*(rand(T-lags-T0,1) > 0.5) - 1;    

            %bootstrapped proxies
            M_boot = M.*(E*ones(1,r));

            %compute the bootstrapped residuals
            U_boot = U_hat(:,:,order).*(E*ones(1,K));

            %produce the bootstrapped left- and right-hand side data
            LHS_boot = zeros(T-lags-T0,K);
            RHS_boot = [ones(T-lags-T0,1),zeros(T-lags-T0,K*lags)];
            RHS_boot(1,:) = RHS(1,:,order);
            for j = 1:T-lags-T0
                LHS_boot(j,:) = RHS_boot(j,:)*A_hat(:,:,order) + U_boot(j,:);
                RHS_boot(j+1,:) = [1,LHS_boot(j,:),RHS_boot(j,2:K*(lags-1)+1)];
            end
            RHS_boot = RHS_boot(1:end-1,:);
        end

        %bootstrap VAR coefficients and innovations
        A_hatboot = (RHS_boot'*RHS_boot)\(RHS_boot'*LHS_boot);
        U_hatboot = LHS_boot - RHS_boot*A_hatboot;

        %identification
        H1_hatboot(:,:,boot,order) = IdentifyFunc(U_hatboot,M_boot);
        
        %generate the IRFs
        IRF_boot(1:K,1,boot,order) =...
            H1_hatboot(:,order,boot,order);
        history_boot = [IRF_boot(1:K,1,boot,order);zeros((lags-1)*K,1)];
        for j = 2:nImp
            IRF_boot(1:K,j,boot,order) =...
                A_hatboot(2:K*lags+1,:)'*history_boot;
            history_boot =...
                [IRF_boot(1:K,j,boot,order);history_boot(1:(lags-1)*K,1)];
        end
        
        %count the number proxy variables not censored to zero
        M_count(boot,:,order) = sum(abs(M_boot) > 0,1);
    end
end

%sort the bootstrapped IRFs
IRF_bootSort = zeros(K,nImp,nBoot,2);
for order = 1
    for n = 1:K
        for j = 1:nImp
            IRF_bootSort(n,j,:,order) = sort(IRF_boot(n,j,:,order));
        end
    end
end

%percentil confidence interval
IRF_confidence = zeros(K,nImp,2,2);
IRF_confidence(:,:,1,1) =...
    IRF_bootSort(:,:,round(nBoot*(1 - alpha)/2),1);
IRF_confidence(:,:,2,1) =...
    IRF_bootSort(:,:,round(nBoot*(1 - (1 - alpha)/2)),1);
IRF_confidence(:,:,1,2) =...
    IRF_bootSort(:,:,round(nBoot*(1 - alpha)/2),2);
IRF_confidence(:,:,2,2) =...
    IRF_bootSort(:,:,round(nBoot*(1 - (1 - alpha)/2)),2);

%display the fewest non-censored proxy observations in the bootstrap
disp([min(M_count(:,:,1));min(M_count(:,:,2))])


%==========================================================================
% Plot IRFs and Confidence Intervals
%==========================================================================
SVAR.LtildeFull = zeros(K,nImp,5,K);
SVAR.LtildeFull(:,:,1,1) = IRF_confidence(:,:,1,1);
SVAR.LtildeFull(:,:,5,1) = IRF_confidence(:,:,2,1);
SVAR.LtildeFull(:,:,3,1) = IRF(:,:,1);

if varSize == 5
    save(strcat('Result5eq_bootstrap','toggleBoot', ...
           num2str(toggleBoot),'p_',num2str(lags)),'SVAR')
    
    fig = figure(1);
    subplot(3,2,1)
    plot(zeros(1,nImp),'k-')
    hold on
    plot(IRF(1,:,1)*100,'b-')
    plot(IRF_confidence(1,:,1,1)*100,'b--')
    plot(IRF_confidence(1,:,2,1)*100,'b--')
    hold off
    axis([0 48 -30 50])
    xlabel('Month')
    title('Federal Funds Rate')

    subplot(3,2,2)
    plot(zeros(1,nImp),'k-')
    hold on
    plot(IRF(2,:,1),'b-')
    plot(IRF_confidence(2,:,1,1),'b--')
    plot(IRF_confidence(2,:,2,1),'b--')
    hold off
    axis([0 48 -1 1])
    xlabel('Month')
    title('Industrial Production')

    subplot(3,2,3)
    plot(zeros(1,nImp),'k-')
    hold on
    plot(IRF(3,:,1)*100,'b-')
    plot(IRF_confidence(3,:,1,1)*100,'b--')
    plot(IRF_confidence(3,:,2,1)*100,'b--')
    hold off
    axis([0 48 -10 15])
    xlabel('Month')
    title('Unemployment Rate')

    subplot(3,2,4)
    plot(zeros(1,nImp),'k-')
    hold on
    plot(IRF(4,:,1),'b-')
    plot(IRF_confidence(4,:,1,1),'b--')
    plot(IRF_confidence(4,:,2,1),'b--')
    hold off
    axis([0 48 -1 0.5])
    xlabel('Month')
    title('PPI')

    subplot(3,2,5)
    plot(zeros(1,nImp),'k-')
    hold on
    plot(IRF(5,:,1)*100,'b-')
    plot(IRF_confidence(5,:,1,1)*100,'b--')
    plot(IRF_confidence(5,:,2,1)*100,'b--')
    hold off
    axis([0 48 -10 15])
    xlabel('Month')
    title('Baa Spread')


    dim = [8,8];
    set(gcf,'paperpositionmode','manual','paperunits','inches');
    set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
    print(fig,'-dpdf',strcat('irf_mon5eq','toggleBoot',num2str(toggleBoot),'p_',num2str(lags)));

elseif varSize == 4
    save(strcat('Result4eq_bootstrap','toggleBoot', ...
                num2str(toggleBoot),'p_',num2str(lags)),'SVAR')
    fig = figure(1);
    subplot(3,2,1)
    plot(zeros(1,nImp),'k-')
    hold on
    plot(IRF(1,:,1)*100,'b-')
    plot(IRF_confidence(1,:,1,1)*100,'b--')
    plot(IRF_confidence(1,:,2,1)*100,'b--')
    hold off
    axis([0 48 -30 50])
    xlabel('Month')
    title('Federal Funds Rate')

    subplot(3,2,2)
    plot(zeros(1,nImp),'k-')
    hold on
    plot(IRF(2,:,1),'b-')
    plot(IRF_confidence(2,:,1,1),'b--')
    plot(IRF_confidence(2,:,2,1),'b--')
    hold off
    axis([0 48 -1 1])
    xlabel('Month')
    title('Industrial Production')

    subplot(3,2,3)
    plot(zeros(1,nImp),'k-')
    hold on
    plot(IRF(3,:,1)*100,'b-')
    plot(IRF_confidence(3,:,1,1)*100,'b--')
    plot(IRF_confidence(3,:,2,1)*100,'b--')
    hold off
    axis([0 48 -10 15])
    xlabel('Month')
    title('Unemployment Rate')

    subplot(3,2,4)
    plot(zeros(1,nImp),'k-')
    hold on
    plot(IRF(4,:,1),'b-')
    plot(IRF_confidence(4,:,1,1),'b--')
    plot(IRF_confidence(4,:,2,1),'b--')
    hold off
    axis([0 48 -1 0.5])
    xlabel('Month')
    title('PPI')

    dim = [8,8];
    set(gcf,'paperpositionmode','manual','paperunits','inches');
    set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
    print(fig,'-dpdf',strcat('irf_mon4eq','toggleBoot',num2str(toggleBoot),'p_',num2str(lags)));

end

toc;

