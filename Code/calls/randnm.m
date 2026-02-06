function store = randnm(med,cov,n)
% Generates pseudo-random numbers with multivariate normal distribution
% with mean vector 'med' and a covariance matrix 'cov':
%
% store (Nxn)   = randnm(med,cov,n)
% med (Nx1)     = mean
% cov (NxN)     = covariance matrix
% n (1x1)       = number of vectors

N = length(cov);
store = zeros(N,n);
for i = 1:n;
    y = randn(N,1);
    store(:,i) = (chol(cov)')*y + med;
end

