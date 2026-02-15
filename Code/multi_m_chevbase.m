%% MASTER-aligned (multi m*): Data processing + VARX-SVAR (Cholesky) + IRFs for multiple m*
%
% Method aligned to master_code.m:
% 1) CMA on full sample
% 2) Cutoff on smoothed data
% 3) Chebyshev basis built on cut sample using chevyp_ort* functions
% 4) VARX + Uhlig-style draws using xx_cond (Schur complement)
% 5) Plot median IRFs for a grid of m*

clear all; clc; close all; tic;
rng(12345,'twister');
addpath("C:\Users\Carlos Coronel\Documents\Informalidad\Informality-InterestRate\Code\calls");

%% 1) Parameters

default_var_indices       = [1,2,3,4];
USE_INTERACTIVE_SELECTION = true;
% 065438
use_precovid   = true;
monthly_data   = true;
USE_UNIT_SHOCK = true;

% Data configuration
if monthly_data
    CUTOFF = "06/01/2019";
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\mensuales_2005.xlsx";
    time_granularity = "Mensual";
else
    CUTOFF = "2019Q2";
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\trimestrales_primero.xlsx";
    time_granularity = "Trimestral";
end

% Grid of m* to compare
m_star_grid = [5 6 7 8 15];

% VAR / IRF
p    = 1;
s    = 49;
reps = 1000;

assert(all(m_star_grid >= 1), 'm_star_grid must contain integers >= 1.');

%% 2) Load and variable selection

TBL = readtable(filename);

% Period formatting (same logic as master_code)
if monthly_data
    try
        if isdatetime(TBL{:,1})
            periods = string(datetime(TBL{:,1}, 'Format', 'MM/dd/yyyy'));
        else
            periods = string(TBL{:,1});
        end
    catch
        warning('Hubo un problema formateando las fechas. Se usaran tal cual vienen.');
        periods = string(TBL{:,1});
    end
else
    periods = string(TBL{:,1});
end

Ytbl_raw    = TBL(:,2:end);
allVarNames = Ytbl_raw.Properties.VariableNames;

if USE_INTERACTIVE_SELECTION
    fprintf('Data Columns:\n');
    for i = 1:numel(allVarNames)
        fprintf('%d: %s\n', i, allVarNames{i});
    end
    sel_idx = input('Indices para ENDOGENAS (VAR / Cholesky Order) (ej: [1, 3]): ');
    fprintf('\nSelecciona las variables EXOGENAS adicionales (ademas del polinomio).\n');
    fprintf('Deja vacio y pulsa Enter si no quieres ninguna.\n');
    sel_exo_idx = input('Indices para EXOGENAS (ej: [5, 6] o []): ');
else
    sel_idx = default_var_indices;
    sel_exo_idx = [];
end

% Endogenous and exogenous selection
Ytbl      = Ytbl_raw(:, sel_idx);
var_order = Ytbl.Properties.VariableNames;

if ~isempty(sel_exo_idx)
    Ytbl_exo_raw = Ytbl_raw(:, sel_exo_idx);
else
    Ytbl_exo_raw = table();
end

if ~isempty(sel_exo_idx)
    external_exo_names = Ytbl_raw.Properties.VariableNames(sel_exo_idx);
else
    external_exo_names = {};
end
if isempty(external_exo_names)
    txt_exo = "Ninguna";
else
    txt_exo = strjoin(string(external_exo_names), ', ');
end

fprintf('VAR variables in final order:\n');
disp(var_order');

[n_full, k_selected] = size(Ytbl);

%% 3) CMA (FULL SAMPLE) -> cutoff -> Chebyshev-ready sample

% 3.1 Full-sample CMA
if monthly_data
    if n_full < 13
        error('Not enough observations for Monthly CMA (Need > 12).');
    end
    idx_cma  = 7:(n_full-6);
    cma_func = @cma_m;
else
    if n_full < 5
        error('Not enough observations for Quarterly CMA (Need > 4).');
    end
    idx_cma  = 3:(n_full-2);
    cma_func = @cma_q;
end

n_cma = numel(idx_cma);

cma_endo = NaN(n_cma, k_selected);
for j = 1:k_selected
    y = double(Ytbl{:,j});
    y_cma = cma_func(y);
    if numel(y_cma) ~= n_cma
        error('CMA output length mismatch in endogenous column %d.', j);
    end
    cma_endo(:,j) = y_cma;
end

cma_exo = [];
if ~isempty(sel_exo_idx)
    n_exo_selected = numel(sel_exo_idx);
    cma_exo = NaN(n_cma, n_exo_selected);
    for j = 1:n_exo_selected
        y = double(Ytbl_exo_raw{:,j});
        y_cma = cma_func(y);
        if numel(y_cma) ~= n_cma
            error('CMA output length mismatch in exogenous column %d.', j);
        end
        cma_exo(:,j) = y_cma;
    end
end

endo_raw    = cma_endo;
periods_cma = periods(idx_cma);

% 3.2 Cutoff after smoothing (same order as master_code)
if use_precovid
    fprintf('\n--- Buscando fecha de corte: %s ---\n', CUTOFF);
    cutoff_idx = find(strcmp(periods_cma, CUTOFF), 1, 'first');

    if isempty(cutoff_idx)
        warning('La fecha "%s" no se encontro. Se usara muestra completa.', CUTOFF);
    else
        endo_raw    = endo_raw(1:cutoff_idx, :);
        periods_cma = periods_cma(1:cutoff_idx);
        if ~isempty(cma_exo)
            cma_exo = cma_exo(1:cutoff_idx, :);
        end
        n_cma = cutoff_idx;
        fprintf('Datos cortados exitosamente en %s. Nueva longitud: %d obs.\n', CUTOFF, cutoff_idx);
    end
end

%% 4) Logs and fixed components

% Endogenous
endo_log = log(endo_raw);

% External exogenous
exo_external_final = [];
if ~isempty(cma_exo)
    if any(cma_exo(:) <= 0)
        warning('Cuidado: Hay valores <= 0 en las exogenas externas antes del log.');
    end
    exo_external_final = log(cma_exo);
end

endo = double(endo_log);

if any(isnan(endo(:))) || any(isnan(exo_external_final(:)))
    error('NaNs detectados. Revisa el CMA o los logs.');
end

%% 5) VARX-Uhlig for multiple m* (MASTER-consistent)

data = endo;
[~, k] = size(data);
ne = k;

% Select one shock index for all m*
variable_names = cellstr(var_order);
fprintf('\nVAR variables en orden final:\n');
for i = 1:numel(variable_names)
    fprintf('%d: %s\n', i, variable_names{i});
end

IMP = 0;
while isempty(IMP) || ~isnumeric(IMP) || IMP < 1 || IMP > k
    try
        IMP = input(sprintf('Elige el indice del shock (1-%d): ', k));
        if isempty(IMP) || ~isnumeric(IMP) || IMP < 1 || IMP > k
            fprintf('Error: debe ser entero entre 1 y %d.\n', k);
        end
    catch
        fprintf('Entrada invalida.\n');
    end
end
fprintf('Shock seleccionado: %d: %s\n', IMP, variable_names{IMP});

nm = numel(m_star_grid);
IRF_med_all = NaN(s, ne, nm);

for idx_m = 1:nm
    m_star = m_star_grid(idx_m);
    fprintf('\n--- Estimando VARX para m* = %d ---\n', m_star);

    % Chebyshev basis exactly as in master_code for this m*
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

    exo = [double(T_cheb), double(exo_external_final)];
    n_exo_total = size(exo, 2);

    fprintf('VARX exogenas totales: %d (= %d Chebyshev + %d externas)\n', ...
        n_exo_total, size(T_cheb,2), size(exo_external_final,2));

    [betas, xxx, omega, res] = estimavarx(data, p, exo);

    betas_imp = betas(:, 1:end-n_exo_total);

    % Conditional X'X for [const + lags] given exogenous (Schur complement)
    r = 1 + k*p;
    q = size(xxx,1);
    if q ~= r + n_exo_total
        error('Dimension mismatch: GAM size does not equal (1+k*p+n_exo_total).');
    end

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

    IRF_draws = NaN(s, ne, reps);

    for rep = 1:reps
        % Draw Sigma
        R = randnm(0, inv(omega)/size(res,1), size(res,1));
        Sigma = inv(R * R');

        % Draw B | Sigma using xx_cond (MASTER-consistent)
        B_cov = kron(inv(xx_cond), Sigma);
        B = randnm(betas_imp(:), B_cov, 1);

        % Rearrange B
        ss_idx = (1:k:length(B)+k)';
        storeb = zeros(k, length(ss_idx)-1);
        for h = 1:length(ss_idx)-1
            storeb(:,h) = B(ss_idx(h):ss_idx(h+1)-1,:);
        end

        % Cholesky shock
        C = transpose(chol(Sigma));
        v_imp = C(:, IMP);
        if USE_UNIT_SHOCK
            scale = v_imp(IMP);
            if abs(scale) > 1e-12
                v_imp = v_imp / scale;
            else
                warning('Shock %d: escala casi cero, se deja como 1 s.d.', IMP);
            end
        end

        resp = impulse(storeb, v_imp, s, p);
        IRF_draws(:,:,rep) = resp;
    end

    IRF_med_all(:,:,idx_m) = median(IRF_draws, 3);
end

%% 6) Plot: multiple m* line IRFs (no confidence bands)

if monthly_data
    freq_divisor = 12;
else
    freq_divisor = 4;
end
x = (0:s-1)/freq_divisor;
years = 0:floor((s-1)/freq_divisor);

colors = lines(nm);

figure('Name','IRFs VARX con Chebyshev (multi m*)','Color','w');

if use_precovid
    cutoff_text = sprintf('Corte en %s', CUTOFF);
else
    cutoff_text = 'Sin corte';
end

sgtitle(sprintf('IRFs (shock: %s) | p=%d, s=%d, %s, %s, Exo=%s', ...
    variable_names{IMP}, p, s, cutoff_text, time_granularity, txt_exo), ...
    'FontSize', 13);

ncols = 4;
nrows = max(2, ceil(ne/ncols));

for ii = 1:(nrows*ncols)
    subplot(nrows, ncols, ii);

    if ii > ne
        axis off;
        continue
    end

    hold on; box on;

    for idx_m = 1:nm
        plot(x, IRF_med_all(:,ii,idx_m), 'LineWidth', 1.4, 'Color', colors(idx_m,:));
    end

    yline(0,'k');
    title(variable_names{ii}, 'Interpreter','none');
    xlabel('Years');
    xticks(years);
    xticklabels(string(years));

    if ii == 1
        legend("m*="+string(m_star_grid), 'Location','best');
    end

    set(gca,'Layer','top');
    hold off;
end

%% 7) Save results

results_meta = struct();
results_meta.var_order           = var_order;
results_meta.exo_names           = external_exo_names;
results_meta.p                   = p;
results_meta.s                   = s;
results_meta.reps                = reps;
results_meta.k                   = k;
results_meta.m_star_grid         = m_star_grid;
results_meta.periods             = periods_cma;
results_meta.shock_plotted_index = IMP;
results_meta.shock_plotted_name  = variable_names{IMP};
results_meta.use_precovid        = use_precovid;
results_meta.cutoff              = CUTOFF;

save('varx_irf_results_multi_mstar_chevbase.mat', 'IRF_med_all', 'results_meta');

toc;
fprintf('Saved as "varx_irf_results_multi_mstar_chevbase.mat".\n');
