function [H1] = IdentifyFunc(U,M)
%==========================================================================
% This function identifies the relevant columns of H using proxy variables.
%   The inputs are:
%       U, the matrix of VAR residuals
%       M, the matrix of proxy variables
%   The outputs are:
%       H1, the columns of interest
%==========================================================================

    %dimensions of the data
    [T,K] = size(U);
    r = size(M,2);

    %covariance matrices
    SigmaU = U'*U/T;
    SigmaU11 = SigmaU(1:r,1:r);
    SigmaU21 = SigmaU(r+1:K,1:r);
    SigmaU22 = SigmaU(r+1:K,r+1:K);

    SigmaMU = M'*U/T;
    SigmaMU1 = SigmaMU(1:r,1:r);
    SigmaMU2 = SigmaMU(1:r,r+1:K);

    %identification
    H21H11 = (SigmaMU1\SigmaMU2)';
    Z = SigmaU22 - H21H11*SigmaU21' -...
        SigmaU21*H21H11' + H21H11*SigmaU11*H21H11';
    H12H12 = (SigmaU21 - H21H11*SigmaU11)'*(Z\(SigmaU21 - H21H11*SigmaU11));
    H11H11 = SigmaU11 - H12H12;
    H22H22 = SigmaU22 - H21H11*H11H11*H21H11';
    H12H22 = (SigmaU21' - H11H11*(H21H11)')/H22H22;

    %estimate B1
    S1S1 = (eye(r) - H12H22*H21H11)*H11H11*(eye(r) - H12H22*H21H11)';
    S1 = chol(S1S1,'lower');
    H11 = (eye(r) - H12H22*H21H11)\S1;
    H1 = [eye(r);H21H11]*H11;
end

