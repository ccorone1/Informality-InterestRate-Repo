function [B,lresm] = svarbq(betas,omega)
% Consider the relation: Ae = Bu: 
% where u is the 'structural' innovation and e is the reduced-form shock
% A is the identity matrix
% [B,lresm] = svarbq(betas,omega)
% 
% betas         = estimated coefficients from the VAR reduced-form
% estimation
% omega (NxN)   = estimated variance-covariance matrix
% B (NxN)       = estimated matrix in the above relation
% lresm (NxN)   = long-run response pattern
% 

betas_c = betas(:,2:size(betas,2));
k = size(betas,1);
p = size(betas,2)/k;

comp1 = eye(k);
for i = 1:p
    comp1 = comp1-betas_c(:,(i-1)*k+1:i*k);
end

B = comp1*chol(comp1\omega/comp1')';
lresm = comp1\B;
