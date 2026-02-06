%**************************************************************************
% ESTIMAVARX.M

% Written by Gustavo Leyva
% Department of Economics
% University of Minnesota
% Created 09.12.2011
% Modified 06.22.2012

% USAGE: [B,GAM,SUUB,U] = estimavarx(data,p,x);
% 
% data : the data set '[y1 y2 y3]'
% p    : the VAR lag order (maximum = 12)
% x    : the set of exogenous variables (like trend variables)
% B    : the matrix of VAR estimated coefficients (by OLS)
% GAM  : matrix X'X
% SUUB : the variance-covariance matrix of the residuals. This matrix is
% degrees of freedom adjuted (dof). If want to work with the unadjusted
% version of that matrix choose SUB
% U    : the vector of estimated residuals
%*************************************************************************/

function [B,GAM,SUUB,U] = estimavarx(data,p,x)

[T,k] = size(data);
T = T-p;
Y = data';
X = x';
cY = size(Y,2);
YY = Y(:,p+1:cY);
[nX,cX] = size(X);
X = X(:,p+1:cX);
Z = zeros( k*p + nX+1 ,T);

for i = 1:T;
    if p == 1;
        Z(:,i) = [1;Y(:,i+p-1);X(:,i)];
    elseif p == 2;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);X(:,i)];
    elseif p == 3;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);X(:,i)];
    elseif p == 4;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);X(:,i)];
    elseif p == 5;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);X(:,i)];
    elseif p == 6;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);X(:,i)];
    elseif p == 7;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);X(:,i)];
    elseif p == 8;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);X(:,i)];
    elseif p == 9;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);X(:,i)];
    elseif p == 10;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);Y(:,i+p-10);X(:,i)];
    elseif p == 11;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);Y(:,i+p-10);Y(:,i+p-11);X(:,i)];
    elseif p == 12;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);Y(:,i+p-10);Y(:,i+p-11);Y(:,i+p-12);X(:,i)];
    end;
end;

B = YY*Z'/(Z*Z');
U = YY-(YY*Z'/(Z*Z'))*Z; U = U';
SUB = (1/T)*YY*(eye(T)-Z'/(Z*Z')*Z)*YY';
SUUB = (T/(T-k*p-1-size(x,2)))*SUB; %bias correction
%SUUB = SUB; %bias correction
GAM = Z*Z';