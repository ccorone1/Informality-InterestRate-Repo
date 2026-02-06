clear; clc;
% reading data

chdir T:\DiscoHDE1\Informality_bs\replication_materials_afterBM\main_text\tables
%AS = xlsread('Table7_compa_ENE','series_empalmadas','B2:S122');
%%
%AS = xlsread('Table7_compa_ENE','series_empalmadas','C15:Q134');
AS = xlsread('Table7_compa_ENE','series_empalmadas','C3:S134');
chdir T:\Carlos_student\LECU
w = 1600;
[~,k] = size(AS);
A = [];
for ii = 1:k
    A(:,ii) = cma_q(AS(:,ii));
end
A = packr(A);
logA = log(A(:,1:k));
A = logA;
T = length(A);
tt=(1:1:T)';
t = 1;
trend = (1:length(A))'; 
trend2 = trend.^2;
trend3 = trend.^3;

su = [];
for j = 1:3
    su = [su chevyp_ort2(j,T,tt) chevyp_ort3(j,T,tt)];
end
%endo = [A(t:T,15) A(t:T,1) A(t:T,3) A(t:T,10) A(t:T,13)  A(t:T,7)]; % igual al orden que aparece en Eviews
endo = [A(t:T,16) A(t:T,1) A(t:T,3) A(t:T,10) A(t:T,2) A(t:T,12)  A(t:T,7)]; % igual al orden que aparece en Eviews
Xt = [chevyp_ort1(T,tt) su];
Y = endo(:,7); %%Endogena a filtrar
X =  [chevyp_ort0(T) chevyp_ort1(T,tt) su];X_T = X';
beta = (X_T*X)\(X_T*Y);
YS = X*beta;
figure(1)
plot([Y YS]);
figure(2)
plot((Y-YS)*100);