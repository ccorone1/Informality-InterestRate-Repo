function [P] = chevyp_ort2(k,n,t)

if k == 1
    X0 = [ones(n,1) (t./n)];X_T0 = X0';
    beta = (X_T0*X0)\(X_T0*chevyp(2*k-1,n,t));
    P = chevyp(2*k-1,n,t) - beta(1) -beta(2)*(t./n);
    P = P/sqrt(mean(P.^2));    
else
    su = [];
    for j = 1:k-1
       su = [su chevyp(2*j-1,n,t)];
    end
    [~,ks] = size(su);
    X0 = [ones(n,1) su (t./n)];X_T0 = X0';
    beta = (X_T0*X0)\(X_T0*chevyp(2*k-1,n,t));
    P = chevyp(2*k-1,n,t) - beta(1) - su*beta(2:1+ks) - beta(2+ks)*(t./n);
    P = P/sqrt(mean(P.^2));
end

end