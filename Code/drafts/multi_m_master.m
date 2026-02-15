
%% Data processing and SVAR with Chebyshev (multiple m*)
%
% Versión: multi-m* (Elizondo)

clear all; clc; close all; tic;
addpath("C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\code\matlab");
addpath("C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\code\matlab\20_11_25");
 

%% 1. Parámetros generales

default_var_indices   = [1,2,3,4];

USE_INTERACTIVE_SELECTION = true;
use_precovid      = true;
use_cosine_def    = true;
monthly_data      = true;
USE_UNIT_SHOCK    = true;

if monthly_data
    CUTOFF   = "12/01/2019";
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\mensuales_2005.xlsx";
    time_granularity = "Mensual";
else
    CUTOFF   = "2019Q4";
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\trimestrales_primero.xlsx";
    time_granularity = "Trimestral";
end

m            = 301;                          % orden máximo Chebyshev
m_star_grid  = [2 3 5 15 50 99 100 149];       % grados a comparar
p            = 1;                            % lags VAR
s            = 84;                           % horizonte IRF
reps         = 1000;                         % draws Uhlig

%% 2. Carga y selección

data = readtable(filename);

% Formato de periodos
if monthly_data
    try
        if isdatetime(data{:,1})
            periods = string(datetime(data{:,1}, 'Format','MM/dd/yyyy'));
        else
            periods = string(data{:,1});
        end
    catch
        warning('Problema formateando fechas mensuales. Se usan como texto.');
        periods = string(data{:,1});
    end
else
    % Columna tipo '2005Q1'
    raw_periods = string(data{:,1});
    try
        dt = datetime(raw_periods,"InputFormat","uuuu'Q'Q");
    catch
        dt = datetime(raw_periods);
    end
    dt = dateshift(dt,'end','quarter');   % fin de trimestre
    periods = string(dt);
end

Ytbl_raw    = data(:,2:end);
allVarNames = Ytbl_raw.Properties.VariableNames;

% Selección interactiva de endógenas / exógenas
if USE_INTERACTIVE_SELECTION
    fprintf('Data Columns:\n');
    for i = 1:numel(allVarNames)
        fprintf('%d: %s\n', i, allVarNames{i});
    end
    sel_idx     = input('Indices para ENDÓGENAS (VAR / Cholesky Order): ');
    sel_exo_idx = input('Indices para EXÓGENAS adicionales ([] si ninguna): ');
else
    sel_idx     = default_var_indices;
    sel_exo_idx = [];
end

Ytbl      = Ytbl_raw(:, sel_idx);
var_order = Ytbl.Properties.VariableNames;

if ~isempty(sel_exo_idx)
    Ytbl_exo_raw = Ytbl_raw(:, sel_exo_idx);
else
    Ytbl_exo_raw = table();
end

% Corte pre-Covid (opcional)
if use_precovid
    cutoff_idx = find(strcmp(periods, CUTOFF));
    if isempty(cutoff_idx)
        warning('Cut off %s no encontrado. Se usa toda la muestra.', CUTOFF);
    else
        Ytbl  = Ytbl(1:cutoff_idx,:);
        if ~isempty(Ytbl_exo_raw)
            Ytbl_exo_raw = Ytbl_exo_raw(1:cutoff_idx,:);
        end
        periods = periods(1:cutoff_idx);
        fprintf('Muestra recortada a %s (obs %d).\n', CUTOFF, cutoff_idx);
    end
end

[n, p_selected] = size(Ytbl);

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

%% 3. CMA y polinomios Chebyshev (base)

if monthly_data
    if n < 13, error('Pocas obs para CMA mensual (requiere >12).'); end
    idx_cma  = 7:(n-6);
    cma_func = @cma_m;
else
    if n < 5, error('Pocas obs para CMA trimestral (requiere >4).'); end
    idx_cma  = 3:(n-2);
    cma_func = @cma_q;
end

n_cma      = numel(idx_cma);
cma_series = NaN(n_cma, p_selected);

for j = 1:p_selected
    y   = double(Ytbl{:,j});
    y_c = cma_func(y);
    if numel(y_c) ~= n_cma
        error('CMA mismatch en col %d.', j);
    end
    cma_series(:,j) = y_c;
end

% CMA para exógenas
if ~isempty(Ytbl_exo_raw)
    n_exo_selected  = length(sel_exo_idx);
    cma_exo_series  = NaN(n_cma, n_exo_selected);
    for j = 1:n_exo_selected
        y   = double(Ytbl_exo_raw{:,j});
        y_c = cma_func(y);
        cma_exo_series(:,j) = y_c;
    end
else
    cma_exo_series = [];
end

endo_raw   = cma_series;
periods_cma = periods(idx_cma);

% Fechas para t_cma
if monthly_data
    % Mensual: ya las guardaste como 'MM/dd/yyyy'
    try
        dt = datetime(periods_cma,"InputFormat","MM/dd/yyyy");
    catch
        dt = datetime(periods_cma); % fallback
    end
else
    % Trimestral: 'periods' ya son fechas tipo '31-Mar-2005',
    % así que NO intentes leerlas como '2005Q1'
    try
        dt = datetime(periods_cma);  % que MATLAB las detecte solo
    catch
        error('No se pudieron convertir los periodos trimestrales a datetime.');
    end
end

time_years = years(dt - dt(1));
t_cma      = 2*(time_years - min(time_years)) / (max(time_years)-min(time_years)) - 1;
t_cma      = max(-1, min(1, t_cma));

% Matriz T de Chebyshev (0..m)
if use_cosine_def
    theta = acos(t_cma);
    kvec  = 0:m;
    T     = cos(theta * kvec);   % n_cma x (m+1)
else
    T       = zeros(n_cma, m+1);
    T(:,1)  = 1;
    if m >= 1, T(:,2) = t_cma; end
    for kk = 2:m
        T(:,kk+1) = 2*t_cma.*T(:,kk) - T(:,kk-1);
    end
end

%% 4. log(Variables) y construcción de endógenas/exógenas externas

endo_log = log(endo_raw);

exo_external_final = [];
if ~isempty(cma_exo_series)
    if any(cma_exo_series(:) <= 0)
        warning('Hay valores <=0 en exógenas antes de log.');
    end
    exo_external_final = log(cma_exo_series);
end

endo = double(endo_log);   % (n_cma x k)

if any(isnan(endo(:)))
    error('NaNs en endógenas tras CMA/log.');
end

%% 5. VARX (Uhlig) para múltiples m*

data = endo;
[~, k] = size(data);
ne     = k;

% Elegir choque (mismo para todos los m*)
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

assert(max(m_star_grid) <= m, 'max(m_star_grid) no puede exceder m.');

% Para la covarianza de B (no depende de exógenas)
xx = momentxx(data, p);

% Contenedor de IRFs medianas
IRF_all = cell(length(m_star_grid), 1);

for idx_m = 1:length(m_star_grid)
    m_star = m_star_grid(idx_m);
    fprintf('\n--- Estimando VARX para m* = %d ---\n', m_star);

    % Polinomio Chebyshev para este m*
    exo_single = T(:, m_star+1);       % columna k = m_star
    exo_single = exo_single - mean(exo_single);

    % Matriz de exógenas completa: [Chebyshev, otras]
    exo = [double(exo_single), double(exo_external_final)];

    % Estimación VARX
    [betas, ~, omega, res] = estimavarx(data, p, exo);

    n_exo_total = size(exo,2);
    betas_imp   = betas(:,1:end-n_exo_total);  % quitamos exógenas

    % IRFs (draws)
    IRF_draws = zeros(s, ne, reps);

    for rep = 1:reps
        % Sigma
        R     = randnm(0, inv(omega)/size(res,1), size(res,1));
        Sigma = inv(R * R');

        % B | Sigma
        B_cov = kron(inv(xx), Sigma);
        B     = randnm(betas_imp(:), B_cov, 1);

        % Reacomodo B en matrices de coef
        ss_idx  = (1:k:length(B)+k)';
        storeb  = zeros(k, length(ss_idx)-1);
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
            else
                warning('Shock %d: escala ~0, se usa 1 s.d.', IMP);
            end
        end

        resp = impulse(storeb, v_imp, s, p);   % s x ne

        IRF_draws(:,:,rep) = resp;
    end

    % Mediana en la dimensión de draws
    IRF_all{idx_m} = median(IRF_draws, 3);
end

%% 6. Gráfica: IRFs para distintos m*

if monthly_data
    freq_divisor = 12;
else
    freq_divisor = 4;
end
x     = (0:s-1)/freq_divisor;
years = 0:floor((s-1)/freq_divisor);

colors = lines(length(m_star_grid));

figure('Name','IRFs VARX con Chebyshev (multi m*)','Color','w');
if use_precovid
    cutoff_text = sprintf('Corte en %s', CUTOFF);
else
    cutoff_text = 'Sin corte';
end

sgtitle(sprintf('IRFs (shock: %s) | p=%d, s=%d, %s, %s, Exo=%s', ...
    variable_names{IMP}, p, s, cutoff_text, time_granularity, txt_exo), ...
    'FontSize', 13);

for ii = 1:ne
    subplot(ceil(ne/2),2,ii);
    hold on; box on;

    for idx_m = 1:length(m_star_grid)
        IRm = IRF_all{idx_m}(:, ii);      % s x 1
        plot(x, IRm, 'LineWidth', 1.4, 'Color', colors(idx_m,:));
    end

    yline(0,'k');
    title(variable_names{ii}, 'Interpreter','none');
    xlabel('Years');
    xticks(years);
    xticklabels(string(years));

    if ii == 1
        legend("m*="+string(m_star_grid),'Location','best');
    end
end

%% 7. Guardar resultados

results_meta = struct();
results_meta.var_order           = var_order;
results_meta.exo_names           = external_exo_names;
results_meta.p                   = p;
results_meta.s                   = s;
results_meta.reps                = reps;
results_meta.k                   = k;
results_meta.m                   = m;
results_meta.m_star_grid         = m_star_grid;
results_meta.periods             = periods_cma;
results_meta.shock_plotted_index = IMP;
results_meta.shock_plotted_name  = variable_names{IMP};

save('varx_irf_results_multi_mstar.mat', ...
     'IRF_all','results_meta');

toc;
fprintf('Saved as "varx_irf_results_multi_mstar.mat".\n');
