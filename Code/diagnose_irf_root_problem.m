%% Diagnose Root Problem In Quarterly IRFs
% Standalone diagnostic to show how unstable draws widen 95% IRFs.
%
% It replicates the key quarterly pipeline:
% 1) CMA on full sample
% 2) Cut sample at cutoff
% 3) Build Chebyshev exogenous terms
% 4) Estimate VARX and draw IRFs
% 5) Compare IRFs using all draws vs only stable draws

clear; clc; close all;

addpath(fullfile(pwd, "calls"));

%% Configuration
filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\trimestrales_primero.xlsx";
cutoff = "2019Q2";          % clean pre-COVID cutoff for quarterly 2x4 CMA
sel_idx = [1,2,3,4];        % endogenous variables (columns in data(:,2:end))
sel_exo_idx = [];           % optional extra external exogenous variables

m_star = 6;                 % Chebyshev degree setup from master
p = 1;                      % VAR lags
s = 17;                     % IRF horizon
reps = 3000;                % Monte Carlo draws
shock_idx = 1;              % which structural shock
resp_idx = 1;               % which response variable to display

stability_tol = 1e-6;
rng(12345, "twister");

%% 1) Load data
tbl = readtable(filename);
periods = string(tbl{:,1});
Ytbl_raw = tbl(:,2:end);
var_names_all = Ytbl_raw.Properties.VariableNames;

Ytbl = Ytbl_raw(:, sel_idx);
if ~isempty(sel_exo_idx)
    Ytbl_exo_raw = Ytbl_raw(:, sel_exo_idx);
else
    Ytbl_exo_raw = table();
end

[n_full, k] = size(Ytbl);
if n_full < 5
    error("Not enough observations for quarterly CMA.");
end

%% 2) Full-sample quarterly CMA
idx_cma = 3:(n_full-2);
n_cma = numel(idx_cma);

cma_endo = NaN(n_cma, k);
for j = 1:k
    y = double(Ytbl{:,j});
    yc = cma_q(y);
    if numel(yc) ~= n_cma
        error("CMA length mismatch in endogenous variable %d.", j);
    end
    cma_endo(:,j) = yc;
end

cma_exo = [];
if ~isempty(sel_exo_idx)
    n_exo_selected = numel(sel_exo_idx);
    cma_exo = NaN(n_cma, n_exo_selected);
    for j = 1:n_exo_selected
        y = double(Ytbl_exo_raw{:,j});
        yc = cma_q(y);
        if numel(yc) ~= n_cma
            error("CMA length mismatch in exogenous variable %d.", j);
        end
        cma_exo(:,j) = yc;
    end
end

endo_raw = cma_endo;
periods_cma = periods(idx_cma);

%% 3) Cut after smoothing
cut_idx = find(strcmp(periods_cma, cutoff), 1, "first");
if isempty(cut_idx)
    error("Cutoff %s not found in periods_cma.", cutoff);
end

endo_raw = endo_raw(1:cut_idx, :);
periods_cma = periods_cma(1:cut_idx);
if ~isempty(cma_exo)
    cma_exo = cma_exo(1:cut_idx, :);
end
n_cma = cut_idx;

%% 4) Build Chebyshev exogenous block (same as master)
tt = (1:n_cma)';
su = [];
Jmax = ceil(m_star/2);
for j = 1:Jmax
    if (2*j - 1) <= m_star
        su = [su, chevyp_ort2(j, n_cma, tt)];
    end
    if (2*j) <= m_star
        su = [su, chevyp_ort3(j, n_cma, tt)];
    end
end
T_cheb = [chevyp_ort1(n_cma, tt), su];
T_cheb = T_cheb - mean(T_cheb, 1);

if any(endo_raw(:) <= 0)
    error("Endogenous CMA has non-positive values before log.");
end
data = log(endo_raw);

exo_external = [];
if ~isempty(cma_exo)
    if any(cma_exo(:) <= 0)
        error("External exogenous CMA has non-positive values before log.");
    end
    exo_external = log(cma_exo);
end

exo = [double(T_cheb), double(exo_external)];
n_exo_total = size(exo,2);

%% 5) Estimate VARX
[betas, xxx, omega, res] = estimavarx(data, p, exo);
betas_imp = betas(:, 1:end-n_exo_total);

r = 1 + k*p;
if n_exo_total > 0
    GAM11 = xxx(1:r, 1:r);
    GAM12 = xxx(1:r, r+1:end);
    GAM22 = xxx(r+1:end, r+1:end);
    if rcond(GAM22) < 1e-12
        GAM22inv = pinv(GAM22);
    else
        GAM22inv = inv(GAM22);
    end
    xx_cond = GAM11 - GAM12 * GAM22inv * GAM12';
else
    xx_cond = xxx;
end

%% 6) Draw IRFs and diagnose roots
rho = NaN(reps,1);
is_stable = false(reps,1);
ir_all = NaN(s, reps);

for i = 1:reps
    R = randnm(0, inv(omega)/size(res,1), size(res,1));
    Sigma = inv(R * R');

    B_cov = kron(inv(xx_cond), Sigma);
    B = randnm(betas_imp(:), B_cov, 1);

    ss = (1:k:length(B)+k)';
    storeb = zeros(k, length(ss)-1);
    for h = 1:length(ss)-1
        storeb(:,h) = B(ss(h):ss(h+1)-1,:);
    end

    % Companion matrix roots for stability test
    B_lags = storeb(:,2:end);
    if p == 1
        G = B_lags;
    else
        ident = [kron(eye(k), eye(p-1)), zeros(k*(p-1), k)];
        G = [B_lags; ident];
    end
    rho(i) = max(abs(eig(G)));
    is_stable(i) = (rho(i) < (1 - stability_tol));

    C = chol(Sigma, "lower");
    v_imp = C(:, shock_idx);
    if abs(v_imp(shock_idx)) > 1e-12
        v_imp = v_imp / v_imp(shock_idx);
    end

    resp = impulse(storeb, v_imp, s, p);
    ir_all(:,i) = resp(:,resp_idx);
end

ir_stable = ir_all(:, is_stable);

if isempty(ir_stable)
    error("No stable draws were found. Try smaller m_star or more data.");
end

Q_all = quantile(ir_all', [0.025, 0.50, 0.975])';
Q_stb = quantile(ir_stable', [0.025, 0.50, 0.975])';

unstable_share = 100 * (1 - mean(is_stable));
fprintf("Obs after CMA+cutoff: %d\n", n_cma);
fprintf("Endogenous vars: %d | Exogenous terms: %d\n", k, n_exo_total);
fprintf("Unstable draws: %.2f%% (%d/%d)\n", unstable_share, sum(~is_stable), reps);
fprintf("Spectral radius quantiles: 90%%=%.4f | 95%%=%.4f | 99%%=%.4f\n", ...
    quantile(rho,0.90), quantile(rho,0.95), quantile(rho,0.99));

%% 7) Plot diagnostic
x = 0:(s-1);
resp_name = var_names_all{sel_idx(resp_idx)};
shock_name = var_names_all{sel_idx(shock_idx)};

figure("Color","w", "Name","Root Diagnostic - Quarterly IRFs");

subplot(1,2,1);
histogram(rho, 40, "FaceColor", [0.30 0.45 0.75], "EdgeColor","none");
hold on;
xline(1.0, "r--", "LineWidth", 1.5);
xline(1.0 - stability_tol, "k:", "LineWidth", 1.2);
title("Spectral Radius of Draws");
xlabel("max |eigenvalue|");
ylabel("Count");
grid on; box on;
hold off;

subplot(1,2,2);
hold on;
patch([x fliplr(x)], [Q_all(:,1)' fliplr(Q_all(:,3)')], [0.80 0.85 0.95], ...
    "EdgeColor","none", "FaceAlpha",0.55);
plot(x, Q_all(:,2), "Color", [0.25 0.35 0.65], "LineWidth", 1.7);

patch([x fliplr(x)], [Q_stb(:,1)' fliplr(Q_stb(:,3)')], [0.92 0.80 0.80], ...
    "EdgeColor","none", "FaceAlpha",0.40);
plot(x, Q_stb(:,2), "Color", [0.70 0.15 0.15], "LineWidth", 1.7);

yline(0, "k-");
title(sprintf("IRF %s <- shock %s", resp_name, shock_name), "Interpreter","none");
xlabel("Horizon (quarters)");
ylabel("Response");
legend("95% all draws","median all draws","95% stable draws","median stable draws", ...
    "Location","best");
grid on; box on;
hold off;

sgtitle(sprintf("Quarterly root diagnostic | cutoff=%s | m*=%d | unstable=%.1f%%", ...
    cutoff, m_star, unstable_share));

