
% PURPOSE: Generates dummy observations for a Minnesota Prior
% -----------------------------------------------------
% USAGE: vm_dummy
% -----------------------------------------------------
% NOTE: requires to run vm_spec.m
% -----------------------------------------------------
 
% tau    : Overall tightness
% d      : Scaling down the variance for the coefficients of a distant lag
% w      : Number of observations used for obtaining the prior for the 
%          covariance matrix of error terms (usually fixed to 1). 
% lambda : Tuning parameter for coefficients for constant.
% mu     : Tuning parameter for the covariance between coefficients*/


%******************************************************** 
% Import data series                                    *
%*******************************************************/
vm_loaddata

nv      = size(YY,2);     %* number of variables */
nobs    = size(YY,1)-T0;  %* number of observations */




%******************************************************** 
% Dummy Observations                                    *
%*******************************************************/
ext_T0  = 0;                           % I don't know why people set up this
                                       % prior like thisni

if MP ==0
    if strcmp(mmodel, '4eq')
        tau	   =   0.5; 
        d	   =   1;
        w	   =   1;
        lambda =   0.5;	 
        mu	   =   0.5;	 
    else
    tau	   =   0.5; 
    d	   =   3;
    w	   =   1;
    lambda =   0.5;	 
    mu	   =   0.5;	 
    end
elseif MP==1
    if strcmp(mmodel, 'm2levNoSp')
        tau	   =   0.5; 
        d	   =   1.5;
        w	   =   1;
        lambda =   1;	 
        mu	   =   1;	 
    else 
    tau	   =   0.2; 
    d	   =   3.5;
    w	   =   1;
    lambda =   0.1;	 
    mu	   =   0.1;	 
    end
end

YY0     =   YY(1:T0,:);  
ybar    =   mean(YY0)';
sbar    =   std(YY0)' ;
premom  =   [ybar sbar];


% Generate matrices with dummy observations
hyp = [tau; d; w; lambda; mu];
[YYdum, XXdum, breakss] = varprior_h(nv,nlags_,nex_,hyp,premom);


% Actual observations

YYact = YY(T0+1:T0+nobs,:);
XXact = zeros(nobs,nv*nlags_);

i = 1;

while (i <= nlags_)
    XXact(:,(i-1)*nv+1:i*nv) = YY(T0-(i-1):T0+nobs-i,:);
    i = i+1;
end

XXact = [XXact ones(nobs,1)];

%******************************************************** 
% Save Results	                                        *
%*******************************************************/

%if uunix == 1                           
%    sname = (strcat('Results/data',char(file_model),'.mat'));
%else
%    sname = (strcat('Results\data',char(file_model),'.mat'));
%end

%save(sname,'YYdum','XXdum','YYact','XXact','YY','T0','nex_','nlags_','nv','ext_T0','pprior')
