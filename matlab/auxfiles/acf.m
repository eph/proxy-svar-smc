function [Rho] = acf(X,J)
%==========================================================================
% This function computes autocorrelations
%   The inputs are:
%       X, the matrix of data
%       J, the number of autocorrelations to produce
%   The outputs are:
%       Rho, matrix of autocorrelations
%==========================================================================

    %dimensions of the data
    [T,K,n] = size(X);
    
    %average 0f the data;
    Xavg = zeros(T,K,n);
    for ii = 1:n
        preXavg = mean(X(:,:,ii));
        Xavg(:,:,ii) = ones(T,1)*preXavg;
    end

    %covariances
    covariance = zeros(J,K,n);
    Rho = zeros(J,K,n);
    for ii = 1:n
        for jj = 1:J
            covariance(jj,:,ii) = diag((X(jj:T,:,ii) - Xavg(jj:T,:,ii))'*...
                (X(1:T-jj+1,:,ii) - Xavg(1:T-jj+1,:,ii)))'/T;
        end
        Rho(:,:,ii) = covariance(:,:,ii)./(ones(J,1)*covariance(1,:,ii));
    end
end

