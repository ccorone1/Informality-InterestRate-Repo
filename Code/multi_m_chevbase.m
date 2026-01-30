%% MASTER (multi m*): Data processing + VARX-SVAR (Cholesky) + IRFs for multiple m*
%
% - Integra la lógica del "Master": CMA en full sample -> corte -> Chebyshev en muestra cortada
% - Para cada m* en m_star_grid:
%     * usa TODA la base Chebyshev k=1..m* (no solo un polinomio)
%     * estima VARX
%     * hace draws tipo Uhlig y guarda la MEDIANA de las IRFs
% - Grafica múltiples m* en la misma figura (SIN intervalos de confianza)
%
% Requiere en el path:
%   - cma_m / cma_q
%   - estimavarx
%   - randnm
%   - impulse
%
% Carlos Coronel - multi m* integrado

clear all; clc; close all; tic;

addpath("C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\code\matlab");
addpath("C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\code\matlab\20_11_25");

%% 1) Parámetros

default_var_indices       = [1,2,3,4];
USE_INTERACTIVE_SELECTION = true;

use_precovid   = true;   % corta en CUTOFF (después del CMA)
use_cosine_def = true;    % Chebyshev via cos(k*acos(t))
monthly_data   = false;
USE_UNIT_SHOCK = true;    % normaliza para que el shock tenga unidad en su propia variable

% --- Datos ---
if monthly_data
    CUTOFF = "12/01/2019";
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\mensuales_2005.xlsx";
    time_granularity = "Mensual";
else
    CUTOFF = "2019Q4";
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\trimestrales_primero.xlsx";
    time_granularity = "Trimestral";
end

% --- Grid de m* a comparar (EDITA AQUÍ) ---
m_star_grid = [ 2 3 4 6 7 8 9 11 12];      % grados a comparar (Chebyshev k=1..m*)
m_max       = max(m_star_grid); % orden máximo para construir T (0..m_max)

% --- VAR / IRF ---
p    = 1;       % rezagos VAR
s    = 49;      % horizonte IRF
reps = 1000;    % draws Uhlig

%% 2) Carga y selección (endo/exo)

TBL = readtable(filename);

% Periodos
if monthly_data
    try
        if isdatetime(TBL{:,1})
            periods = string(datetime(TBL{:,1}, 'Format', 'MM/dd/yyyy'));
        else
            periods = string(TBL{:,1});
        end
    catch
        warning('Problema formateando fechas mensuales. Se usan tal cual.');
        periods = string(TBL{:,1});
    end
else
    periods = string(TBL{:,1});   % ej: '2005Q1'
end

Ytbl_raw    = TBL(:,2:end);
allVarNames = Ytbl_raw.Properties.VariableNames;

% Selección interactiva
if USE_INTERACTIVE_SELECTION
    fprintf('Data Columns:\n');
    for i = 1:numel(allVarNames)
        fprintf('%d: %s\n', i, allVarNames{i});
    end
    sel_idx = input('Indices para ENDÓGENAS (VAR / Cholesky Order) (ej: [1 3]): ');
    fprintf('\nSelecciona variables EXÓGENAS adicionales (además del polinomio). [] si ninguna.\n');
    sel_exo_idx = input('Indices para EXÓGENAS (ej: [5 6] o []): ');
else
    sel_idx     = default_var_indices;
    sel_exo_idx = [];
end

% Endógenas
Ytbl      = Ytbl_raw(:, sel_idx);
var_order = Ytbl.Properties.VariableNames;

% Exógenas externas (raw)
if ~isempty(sel_exo_idx)
    Ytbl_exo_raw = Ytbl_raw(:, sel_exo_idx);
else
    Ytbl_exo_raw = table();
end

% Nombres exógenas externas
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

%% 3) CMA (FULL SAMPLE) -> corte -> Chebyshev (en muestra final)

% 3.1 CMA sobre muestra completa
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

% CMA endógenas
cma_endo = NaN(n_cma, k_selected);
for j = 1:k_selected
    y = double(Ytbl{:,j});
    y_cma = cma_func(y);
    if numel(y_cma) ~= n_cma
        error('CMA output length mismatch in endo col %d.', j);
    end
    cma_endo(:,j) = y_cma;
end

% CMA exógenas externas
cma_exo = [];
if ~isempty(sel_exo_idx)
    n_exo_selected = numel(sel_exo_idx);
    cma_exo = NaN(n_cma, n_exo_selected);
    for j = 1:n_exo_selected
        y = double(Ytbl_exo_raw{:,j});
        y_cma = cma_func(y);
        if numel(y_cma) ~= n_cma
            error('CMA output length mismatch in exo col %d.', j);
        end
        cma_exo(:,j) = y_cma;
    end
end

endo_raw    = cma_endo;
periods_cma = periods(idx_cma);

% 3.2 Corte (ya suavizado)
if use_precovid
    fprintf('\n--- Buscando fecha de corte: %s ---\n', CUTOFF);
    cutoff_idx = find(strcmp(periods_cma, CUTOFF));
    if isempty(cutoff_idx)
        warning('La fecha "%s" no se encontró. Se usará muestra completa.', CUTOFF);
    else
        endo_raw    = endo_raw(1:cutoff_idx, :);
        periods_cma = periods_cma(1:cutoff_idx);
        if ~isempty(cma_exo)
            cma_exo = cma_exo(1:cutoff_idx, :);
        end
        n_cma = cutoff_idx;
        fprintf('Datos cortados en %s. Nueva longitud: %d obs.\n', CUTOFF, cutoff_idx);
    end
end

% 3.3 Construcción t in [-1,1] y base Chebyshev 0..m_max
if monthly_data
    try
        dt = datetime(periods_cma, "InputFormat","MM/dd/yyyy");
    catch
        dt = datetime(periods_cma);
    end
else
    try
        dt = datetime(periods_cma, "InputFormat","uuuu'Q'Q");
    catch
        dt = datetime(periods_cma);
    end
    dt = dateshift(dt, 'end', 'quarter');
end

time_years = years(dt - dt(1));
t_cma = 2*(time_years - min(time_years)) / (max(time_years) - min(time_years)) - 1;
t_cma = max(-1, min(1, t_cma));

if use_cosine_def
    theta = acos(t_cma);
    kvec  = 0:m_max;
    T     = cos(theta * kvec);   % n_cma x (m_max+1)
else
    T       = zeros(n_cma, m_max+1);
    T(:,1)  = 1;
    if m_max >= 1, T(:,2) = t_cma; end
    for kk = 2:m_max
        T(:,kk+1) = 2*t_cma.*T(:,kk) - T(:,kk-1);
    end
end

%% 4) Logs y construcción final endo/exo externas (se mantendrán fijas en el loop)

% Endógenas: log(CMA)
if any(endo_raw(:) <= 0)
    error('Hay valores <= 0 en endógenas después del CMA; no se puede aplicar log.');
end
endo = log(endo_raw);
endo = double(endo);

% Exógenas externas: log(CMA_exo) si existen
exo_external_final = [];
if ~isempty(cma_exo)
    if any(cma_exo(:) <= 0)
        warning('Hay valores <= 0 en exógenas externas antes del log.');
    end
    exo_external_final = log(cma_exo);
    exo_external_final = double(exo_external_final);
end

if any(isnan(endo(:))) || any(isnan(exo_external_final(:)))
    error('NaNs detectados tras CMA/log.');
end

%% 5) VARX-Uhlig para múltiples m* (usa base Chebyshev completa 1..m*)

data = endo;
[~, k] = size(data);
ne = k;

% Elegir shock (mismo para todos los m*)
variable_names = cellstr(var_order);
fprintf('\nVar variables en orden final:\n');
for i = 1:numel(variable_names)
    fprintf('%d: %s\n', i, variable_names{i});
end

IMP = 0;
while isempty(IMP) || ~isnumeric(IMP) || IMP < 1 || IMP > k
    try
        IMP = input(sprintf('Elige el índice del shock (1-%d): ', k));
        if isempty(IMP) || ~isnumeric(IMP) || IMP < 1 || IMP > k
            fprintf('Error: debe ser entero entre 1 y %d.\n', k);
        end
    catch
        fprintf('Entrada inválida.\n');
    end
end
fprintf('Shock seleccionado: %d: %s\n', IMP, variable_names{IMP});

assert(max(m_star_grid) <= m_max, 'max(m_star_grid) no puede exceder m_max.');

% Contenedor: mediana IRF por m*  -> (s x ne x nm)
nm = numel(m_star_grid);
IRF_med_all = NaN(s, ne, nm);

for idx_m = 1:nm
    m_star = m_star_grid(idx_m);
    fprintf('\n--- Estimando VARX para m* = %d (Chebyshev k=1..%d) ---\n', m_star, m_star);

    % Base Chebyshev k=1..m_star (drop T0 porque estimavarx mete intercepto)
    T_cheb = T(:, 2:(m_star+1));      % n_cma x m_star
    T_cheb = T_cheb - mean(T_cheb, 1);

    % Exógenas finales para este m*
    exo = [double(T_cheb), double(exo_external_final)];
    n_exo_total = size(exo, 2);

    fprintf('VARX exógenas totales: %d (= %d Chebyshev + %d externas)\n', ...
        n_exo_total, size(T_cheb,2), size(exo_external_final,2));

    % Estimar VARX
    [betas, xxx, omega, res] = estimavarx(data, p, exo);

    % Quitar exógenas de betas (para dinámica de la VAR)
    betas_imp = betas(:, 1:end-n_exo_total);

    % X'X condicional en exógenas (Schur complement) como en el master
    r = 1 + k*p;           % const + lags
    q = size(xxx,1);       % total regresores = r + n_exo_total
    if q ~= r + n_exo_total
        error('Dimension mismatch: size(xxx,1) != (1+k*p+n_exo_total).');
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

    % Draws y mediana (SIN intervalos)
    IRF_draws = NaN(s, ne, reps);

    for rep = 1:reps
        % Sigma ~ invWishart-ish (mismo esquema que tu master)
        R     = randnm(0, inv(omega)/size(res,1), size(res,1));
        Sigma = inv(R * R');

        % B | Sigma
        B_cov = kron(inv(xx_cond), Sigma);
        B     = randnm(betas_imp(:), B_cov, 1);

        % Reacomodo de B -> storeb (k x r)
        ss_idx = (1:k:length(B)+k)';     % r+1 cortes
        storeb = zeros(k, length(ss_idx)-1);
        for h = 1:length(ss_idx)-1
            storeb(:,h) = B(ss_idx(h):ss_idx(h+1)-1,:);
        end

        % Cholesky e impulso
        C     = transpose(chol(Sigma));
        v_imp = C(:, IMP);

        if USE_UNIT_SHOCK
            scale = v_imp(IMP);
            if abs(scale) > 1e-12
                v_imp = v_imp / scale;
            end
        end

        % IRF: s x ne
        resp = impulse(storeb, v_imp, s, p);
        IRF_draws(:,:,rep) = resp;
    end

    IRF_med_all(:,:,idx_m) = median(IRF_draws, 3);
end

%% 6) Gráfica: múltiples m* (sin intervalos)

if monthly_data
    freq_divisor = 12;
    xlabel_text  = 'Years';
else
    freq_divisor = 4;
    xlabel_text  = 'Years';
end
x     = (0:s-1)/freq_divisor;
years = 0:floor((s-1)/freq_divisor);

colors = lines(nm);

figure('Name','IRFs VARX Chebyshev (multi m*)','Color','w');

if use_precovid
    cutoff_text = sprintf('Corte en %s', CUTOFF);
else
    cutoff_text = 'Sin corte';
end

sgtitle(sprintf('IRFs multi m* (shock: %s) | p=%d, s=%d, %s, %s, Exo=%s', ...
    variable_names{IMP}, p, s, cutoff_text, time_granularity, txt_exo), 'FontSize', 13);

% Layout tipo master (2x4) con slots vacíos si ne<8
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
    xlabel(xlabel_text);
    xticks(years);
    xticklabels(string(years));

    if ii == 1
        legend("m*="+string(m_star_grid), 'Location','best');
    end

    set(gca,'Layer','top');
    hold off;
end

%% 7) Guardar resultados

results_meta = struct();
results_meta.var_order           = var_order;
results_meta.exo_names           = external_exo_names;
results_meta.p                   = p;
results_meta.s                   = s;
results_meta.reps                = reps;
results_meta.k                   = k;
results_meta.m_star_grid         = m_star_grid;
results_meta.m_max               = m_max;
results_meta.periods             = periods_cma;
results_meta.shock_plotted_index = IMP;
results_meta.shock_plotted_name  = variable_names{IMP};
results_meta.use_precovid        = use_precovid;
results_meta.cutoff              = CUTOFF;

save('varx_irf_results_multi_mstar_MASTER.mat', 'IRF_med_all', 'results_meta');

toc;
fprintf('Saved as "varx_irf_results_multi_mstar_MASTER.mat".\n');
