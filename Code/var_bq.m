%
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% Contact: gustavo.leyva.jimenez@gmail.com
% Created on June, 2009
% Modified on July 13th 2012
% 
% This program performs estimation and inference on a VAR system using the
% Blanchard-Quah decomposition.
% Parameters are estimated by OLS.
% Inference is based on Uhlig(1998) and Uhlig(2005) using Bayesian
% techniques. This program uses the weak prior following
% Uhlig(2005).
%
% References:
% Uhlig(1998): "The robustness of identified VAR conclusions about money. A
% comment." Carnegie-Rochester Conference Series on Public Policy 49,
% 245-263. North-Holland.
% Uhlig(2005): "What are the effects of monetary policy on output? Results
% from an agnostic identification procedure." Journal of Monetary Economics
% 52, 381-419.
%
tic;
%% reading data (G:\DAT\GIE\SVARs_TCR\Versión_2012\New_Data\base_version2012.xls)
%data = xlsread('base_version2012',2,'A1:H314');
%data = xlsread('base_version2012',5,'A1:P105');         % version data trimestralizada

run transformacion_datos.m

data = endo;                 % ---------------- data finally used                     
p = 1;                       % ---------------- VAR lag order
s = 110;                      % ---------------- impulse-response horizon
k = size(data,2);            % ---------------- number of variables
reps = 100;                 % ---------------- number of Monte Carlo replications
ne = size(endo,2);           % ---------------- number of endo variables
t=(0:s-1)';

[betas,xxx,omega,res] = estimavarx(data,p,exo);         % ----------- OLS estimation                          
betas_imp = betas(:,1:size(betas,2)-size(exo,2));       % ----------- betas without trend's coefficients%
%betas_imp = betas(:,1:size(betas,2));                  % ----------- when not using trend
%-------------------------------------------------------
%% format of 'betas'
%-------------------------------------------------------
% y1 = cons1 a1*y1(-1) a2*y2(-1) a3*y3(-1) a4*y4(-1) ...
% y2 = cons2 b1*y1(-1) b2*y2(-1) b3*y3(-1) b4*y4(-1) ...
% y3 = cons3 c1*y1(-1) c2*y2(-1) c3*y3(-1) c4*y4(-1) ...
% y4 = cons4 d1*y1(-1) d2*y2(-1) d3*y3(-1) d4*y4(-1) ...
%-------------------------------------------------------

%% allocating space for the impulse-response functions
for ii = 1:ne
    for rr = 1:ne
        eval(['resp' num2str(rr) 'imp' num2str(ii) '=' 'zeros(s,reps);'  ]);
    end    
end    
%% allocating space for the variance decomposition functions
%var_resp1imp1 = zeros(s,reps);
%var_resp2imp1 = zeros(s,reps);
%var_resp3imp1 = zeros(s,reps);
%var_resp1imp2 = zeros(s,reps);
%var_resp2imp2 = zeros(s,reps);
%var_resp3imp2 = zeros(s,reps);
%var_resp1imp3 = zeros(s,reps);
%var_resp2imp3 = zeros(s,reps);
%var_resp3imp3 = zeros(s,reps);

i = 1;
while i <= reps
    
    %% starting the drawing from the posterior distribution - see Uhlig(1998)    
    R = randnm(0,inv(omega)/size(res,1),size(res,1));
    Sigma = inv(R*R');
    xx = momentxx(data,p);
    B_cov = kron(inv(xx),Sigma);
    B = randnm((betas_imp(:)),B_cov,1);
    %% reshaping B in the format above (result = storeb)
    ss = 1:k:length(B)+k; ss = ss';
    storeb = zeros(k,length(ss)-1);
    for h = 1:length(ss)-1
        storeb(:,h) = B(ss(h):ss(h+1)-1,:);
    end
    %% calling BQ svar procedure
    [b,iaab] = svarbq(storeb,Sigma);
    %% defining the impulse vectors
    for ii = 1:ne
        eval(['imp' num2str(ii) '=' 'b(:,' num2str(ii) ');' ]);        
    end     
    %% constructing the response functions   
    for ii = 1:ne
        eval(['resp' num2str(ii) '=' 'impulse(storeb,imp' num2str(ii) ',s,p);' ]);        
    end
    for ii = 1:ne
        eval(['resp' num2str(ii) '=' 'cumsum(resp' num2str(ii) ');' ]);        
    end    
    %resp1 = cumsum(resp1);            % accumulated responses to shock 1
    %resp2 = cumsum(resp2);            % accumulated responses to shock 2
    %resp3 = cumsum(resp3);            % accumulated responses to shock 3
    %% saving the response/variance decomposition functions
    for ii = 1:ne
        for rr = 1:ne
            eval(['resp' num2str(rr) 'imp' num2str(ii) '(:,' num2str(i) ')' '=' 'resp' num2str(ii) '(:,' num2str(rr) ');' ]);
        end    
    end      
    
    for ii = 1:ne
        for rr = 1:ne
            aux = zeros(s,1);
            for rrr = 1:ne
                eval(['aux = aux + resp' num2str(rrr) '(:,' num2str(rr) ').^2;' ]);
            end
            eval(['var_resp' num2str(rr) 'imp' num2str(ii) '(:,' num2str(i) ')' '=' 'resp' num2str(ii) '(:,' num2str(rr) ').^2./aux;' ]);
        end    
    end          
    i = i+1;
    
end

%% calculating the elasticity
color0=[0.678, 0.847, 0.902];
color1=[0.5843,0.8157,0.9882];
color2=[0,0,0.5];
m=1;
x=t';
IMP=1;
%%
eval(['IR = quantile(transpose(resp1imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(resp1imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(1)
subplot(3,3,1)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(resp2imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(resp2imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(1)
subplot(3,3,2)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(resp3imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(resp3imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(1)
subplot(3,3,3)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(resp4imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(resp4imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(1)
subplot(3,3,4)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(resp5imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(resp5imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(1)
subplot(3,3,5)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(resp6imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(resp6imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(1)
subplot(3,3,6)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(resp7imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(resp7imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(1)
subplot(3,3,7)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
set(gca, 'Layer', 'top')
%plot(t,resp6imp1*100,'LineWidth',0.25);
%hold on
%plot(t,IR2*100,'LineWidth',10);
%hold on
%plot(t,zeros(s,1),'LineWidth',2);
%XX = [t IRm*100];
%YY = [t(2:length(t)) IR(2:length(t),:)*100];
%fanplot(XX,YY,'NumQuantiles',3)
toc
%%
eval(['IR = quantile(transpose(var_resp1imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(var_resp1imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(2)
subplot(3,3,1)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
ylim([0 100])
yticks([0 25 50 75 100])
yticklabels({'0','25','50','75','100'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(var_resp2imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(var_resp2imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(2)
subplot(3,3,2)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
ylim([0 100])
yticks([0 25 50 75 100])
yticklabels({'0','25','50','75','100'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(var_resp3imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(var_resp3imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(2)
subplot(3,3,3)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
ylim([0 100])
yticks([0 25 50 75 100])
yticklabels({'0','25','50','75','100'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(var_resp4imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(var_resp4imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(2)
subplot(3,3,4)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
ylim([0 100])
yticks([0 25 50 75 100])
yticklabels({'0','25','50','75','100'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(var_resp5imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(var_resp5imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(2)
subplot(3,3,5)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
ylim([0 100])
yticks([0 25 50 75 100])
yticklabels({'0','25','50','75','100'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(var_resp6imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(var_resp6imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(2)
subplot(3,3,6)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
ylim([0 100])
yticks([0 25 50 75 100])
yticklabels({'0','25','50','75','100'})
set(gca, 'Layer', 'top')
%
eval(['IR = quantile(transpose(var_resp7imp' num2str(IMP) '),[0.05 0.5 0.95]);IR=transpose(IR);IR=IR*100;' ]);
eval(['IRm = median(var_resp7imp' num2str(IMP) '*100,2);' ]);
y1=IR(:,1)';
y2=IR(:,2)';
y3=IR(:,3)';
figure(2)
subplot(3,3,7)
hold all
plot(x,y1);
plot(x,y2);
plot(x,y3);
patch([x fliplr(x)], [y1 fliplr(y2)], color1,'EdgeColor',color1,'FaceColor',color1)
hold on
patch([x fliplr(x)], [y2 fliplr(y3)], m*color1,'EdgeColor',color1,'FaceColor',m*color1)
%alpha(0.5)
box on
hold on
plot(t,zeros(s,1),'r','LineWidth',1.5);
hold on
plot(x,IRm','color',color2,'LineWidth',2);
box on
xticks([0 2 4 6 8 10])
xticklabels({'0','2','4','6','8','10'})
ylim([0 100])
yticks([0 25 50 75 100])
yticklabels({'0','25','50','75','100'})
set(gca, 'Layer', 'top')
%plot(t,resp6imp1*100,'LineWidth',0.25);
%hold on
%plot(t,IR2*100,'LineWidth',10);
%hold on
%plot(t,zeros(s,1),'LineWidth',2);
%XX = [t IRm*100];
%YY = [t(2:length(t)) IR(2:length(t),:)*100];
%fanplot(XX,YY,'NumQuantiles',3)
toc