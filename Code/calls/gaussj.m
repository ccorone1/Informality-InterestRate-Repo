function x = gaussj(A,b)
%%
% Written by Yousef Saad
%---------------------------------------------------
% function x = gaussj (A,b)
% solves A x = b by Gauss-Jordan elimination
%---------------------------------------------------
n = size(A,1);
A = [A,b];

for k = 1:n
    for i = 1:n
        if (i ~= k)
            piv = A(i,k) / A(k,k);
            A(i,k+1:n+1) = A(i,k+1:n+1) - piv*A(k,k+1:n+1);
        end
    end
end

x = A(:,n+1) ./ diag(A);