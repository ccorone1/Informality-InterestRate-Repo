%% PS1_2
% Written by Gustavo Leyva
clc; clear;

cd('C:\Users\Carlos Coronel\Documents\GitHub\Informality\Informality-InterestRate\Data');

igae = readmatrix('mensuales_2005.xlsx', 'Sheet', 'Hoja1', 'Range', 'B2:B249');

%%
T = length(igae);
tt=(1:1:T)';
su = [];
for j = 1:3
    su = [su chevyp_ort2(j,T,tt) chevyp_ort3(j,T,tt)]; % order seven
end
Y = log(igae);
X =  [chevyp_ort0(T) chevyp_ort1(T,tt) su];
A = (X'*X);
b = (X'*Y);
beta_ols = gaussj(A,b);
YS = X*beta_ols;
figure(1)
plot([Y YS]);
figure(2)
plot((Y-YS)*100);