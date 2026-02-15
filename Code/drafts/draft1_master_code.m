%% Data processing and SVAR with Chevysheb estimation with input selection
%
% By Carlos Coronel (Last update: 27 November, 2025)                                                  
%
% Sections:
%
% 1. Loading, selection and cutoff of data.
% 2. Deseasonalized (CMA) y chebyshev pol.
% 3. Log variables
% 4. Varx & IRF (shock selection)
%
% There are calls to functions created by Gustavo Leyva not included in
% this script, such as:
% - cma_q.m
% - cma_m.m (by Carlos Coronel)
% - estimavarx.m
% - randnm.m
% - momentxx.m
% - impulse.m
%
% This code is thought to work with the following automatizations: you can
% choose if you want to use an interactive selection of the VAR and exo
% variables (for which a Chevyshev Polynomial is used in addition to the
% Time series exogenous variables), you can choose between using cosine or
% recursive definition of the polynomial [results don't change], and choose
% the degree of the polynomial to be used (m_star). You can choose to use 68% or 95%
% CI. And whether to cut or not at 2019Q4 (this is hardcoded).


clear all; clc; close all; tic;
addpath("C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\code\matlab");
addpath("C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\code\matlab\20_11_25");
 

%%  1. Parameters

% default index
default_var_indices = [1,2,3,4];

USE_INTERACTIVE_SELECTION = true; %True: manual selection, false: default

use_precovid = true; %crop to cutoff

use_cosine_def = true; % true: T_k(x)=cos(k*acos(x)); false: recurrencia T_k

use_95percent = false; %true: CI to 95%, false: 68%

monthly_data = true; %Usar datos mensuales

USE_UNIT_SHOCK = true;


% --- Parámetros del Modelo ---

if monthly_data
    CUTOFF = "12/01/2020"; 
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\mensuales_2005.xlsx";
    time_granularity = "Mensual";
else
    CUTOFF = "2020Q4";
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\trimestrales_primero.xlsx";
    time_granularity = "Trimestral";
end
m = 16;      % Chebyshev max order
m_star = 15; % Degree of the polynomial to be used
p = 1;      % VAR lags
s = 84;     % IRF horizon (steps)
reps = 1000;  % Monte carlo repetitions

%% 2. Loading, selection and processing
%  loading
% Usar lectura por defecto de Matlab como solicitado
data = readtable(filename);

% Formatting periods
if monthly_data
    try
        if isdatetime(data{:,1})
             periods = string(datetime(data{:,1}, 'Format', 'MM/dd/yyyy'));
        else
             % Si Matlab lo leyó como texto/string
             periods = string(data{:,1});
        end
    catch
        warning('Hubo un problema formateando las fechas. Se usarán tal cual vienen.');
        periods = string(data{:,1});
    end
else
    periods = string(data{:,1});    % Column A (periods)
end

Ytbl_raw = data(:,2:end);       % Column B (series)
allVarNames = Ytbl_raw.Properties.VariableNames;

% --- Dynamic selection (Endógenas) ---
if USE_INTERACTIVE_SELECTION 
    fprintf('Data Columns:\n');
    for i = 1:numel(allVarNames)
        fprintf('%d: %s\n', i, allVarNames{i});
    end
    sel_idx = input('Indices para ENDÓGENAS (VAR / Cholesky Order) (ej: [1, 3]): ');
    
    % --- NUEVO: Selección de Exógenas ---
    fprintf('\nSelecciona las variables EXÓGENAS adicionales (además del polinomio).\n');
    fprintf('Deja vacío y pulsa Enter si no quieres ninguna.\n');
    sel_exo_idx = input('Indices para EXÓGENAS (ej: [5, 6] o []): ');
else
    sel_idx = default_var_indices;
    sel_exo_idx = []; % Por defecto ninguna si no es interactivo
end

% Filtrar Endógenas
Ytbl = Ytbl_raw(:, sel_idx);
var_order = Ytbl.Properties.VariableNames;

% Filtrar Exógenas (Raw)
if ~isempty(sel_exo_idx)
    Ytbl_exo_raw = Ytbl_raw(:, sel_exo_idx);
else
    Ytbl_exo_raw = table(); % Tabla vacía
end

fprintf('Var variables in final order:\n');
disp(var_order');

% Covid Cutoff
if use_precovid
    cutoff_idx = find(strcmp(periods, CUTOFF));
    if isempty(cutoff_idx)
        warning('Cut off date was not found "%s". All the data will be used.', CUTOFF);
    else
        Tcut = 1:cutoff_idx;
        Ytbl = Ytbl(Tcut, :);
        
        if ~isempty(Ytbl_exo_raw)
            Ytbl_exo_raw = Ytbl_exo_raw(Tcut, :);
        end
        
        periods = periods(Tcut);
        fprintf('Cropped data to %s (line %d).\n', CUTOFF, cutoff_idx);
    end
end
[n, p_selected] = size(Ytbl); % n = obs, p_selected = # de variables elegidas

% Extraer nombres de las exógenas externas ---
if ~isempty(sel_exo_idx)
    % Extrae los nombres de las columnas seleccionadas
    external_exo_names = Ytbl_raw.Properties.VariableNames(sel_exo_idx);
else
    external_exo_names = {};
end

if isempty(external_exo_names)
    txt_exo = "Ninguna";
else
    % Une los nombres con una coma y espacio. Ej: "Tasa, Petroleo"
    txt_exo = strjoin(string(external_exo_names), ', '); 
end
%% 3. CMA & CHEBYSHEV
% CMA Logic Adjustment for Frequency
if monthly_data
    % Monthly: Centered 12-term MA loses 6 obs on each side (indices 7 to n-6)
    if n < 13
        error('Error: Not enough observations for Monthly CMA (Need > 12).');
    end
    idx_cma = 7:(n-6);
    cma_func = @cma_m; % Handle to monthly function
else
    % Quarterly: Centered 4-term MA loses 2 obs on each side (indices 3 to n-2)
    if n < 5
        error('Error: Not enough observations for Quarterly CMA (Need > 4).');
    end
    idx_cma = 3:(n-2);
    cma_func = @cma_q; % Handle to quarterly function
end

n_cma = numel(idx_cma);
cma_series = NaN(n_cma, p_selected);

% Execute CMA for endo
for j = 1:p_selected
    y = double(Ytbl{:,j});
    try
        y_cma = cma_func(y); 
        % Safety check on dimensions
        if length(y_cma) ~= n_cma
            error('CMA output length mismatch. Expected %d, got %d', n_cma, length(y_cma));
        end
        cma_series(:,j) = y_cma;
    catch ME
        warning('Error calculating CMA for column %s: %s', var_order{j}, ME.message);
    end
end


n_exo_selected = length(sel_exo_idx);
cma_exo_series = []; % Inicializar vacío

% Execute CMA for exo
if n_exo_selected > 0
    cma_exo_series = NaN(n_cma, n_exo_selected);
    for j = 1:n_exo_selected
        y = double(Ytbl_exo_raw{:,j}); % Usar la tabla de exógenas
        try
            y_cma = cma_func(y); % Misma función (cma_q o cma_m)
            cma_exo_series(:,j) = y_cma;
        catch ME
            warning('Error en CMA exógena %d: %s', j, ME.message);
        end
    end
end

% cma_exo_series is of the same lenght (n_cma) of that of the endogenous
% variables


% 'endo_raw' deseasonalized series
endo_raw = cma_series;
periods_cma = periods(idx_cma);
% Chevysheb
if monthly_data
    % Fechas mensuales (formato exportado en tus archivos)
    try
        dt = datetime(periods_cma, "InputFormat","MM/dd/yyyy");
    catch
        dt = datetime(periods_cma); % fallback
    end
else
    % Fechas trimestrales tipo "2005Q1"
    try
        dt = datetime(periods_cma, "InputFormat","uuuu'Q'Q");
    catch
        dt = datetime(periods_cma);
    end
    % Ajustar al final del trimestre (económicamente estándar)
    dt = dateshift(dt, 'end', 'quarter');
end

% --- Construir tiempo en años desde inicio ---
time_years = years(dt - dt(1));  % ej. [0, 1/12, 2/12, ..., 5]

% --- Escalar a [-1, 1] ---
t_cma = 2*(time_years - min(time_years)) / (max(time_years) - min(time_years)) - 1;
t_cma = max(-1, min(1, t_cma));  % seguridad

% --- Construir el polinomio ---
if use_cosine_def
    theta = acos(t_cma);
    kvec  = 0:m;
    T     = cos(theta * kvec);
else
    T        = zeros(n_cma, m+1);
    T(:, 1)  = 1;
    if m >= 1
        T(:, 2) = t_cma;
    end
    for kk = 2:m
        T(:, kk+1) = 2*t_cma.*T(:, kk) - T(:, kk-1);
    end
end

% Selección del polinomio
assert(m_star >= 1 && m_star <= m, 'm_star must be between 1 and m.');
exo_single = T(:, m_star+1);
exo_single = exo_single - mean(exo_single);  % centrar

%% 4. log(Variable) y Construcción Final
% Log a las endógenas 
endo_log = log(endo_raw); 

% Log a las exógenas externas (si existen)
exo_external_final = [];
if ~isempty(cma_exo_series)
    % Verificación de positivos
    if any(cma_exo_series(:) <= 0)
        warning('Cuidado: Hay valores <= 0 en las exógenas externas antes del log.');
    end
    exo_external_final = log(cma_exo_series);
end

% Input variables final
endo = double(endo_log);
% CONCATENACIÓN: [Polinomio, Series_Externas] (El orden no altera el VAR)
exo = [double(exo_single), double(exo_external_final)];

% Chequeos
assert(size(endo,1) == size(exo,1), 'Exo y Endo tienen distinto largo temporal.');
if any(isnan(endo(:))) || any(isnan(exo(:)))
    error('NaNs detectados. Revisa el CMA o los logs.');
end

%% 5. VARX (UHLIG)
data = endo; 
[~, k] = size(data); 
ne = k;
t = (0:s-1)';

% assert(size(exo,2)==1, 'UNA sola exógena (columna).'); <-- BORRADO O COMENTADO
fprintf('Estimando VARX con %d variables exógenas (1 Chebyshev + %d Externas)\n', size(exo,2), size(exo_external_final,2));

% VARX con OLS (La función estimavarx debe soportar matriz en 'exo')
[betas, xxx, omega, res] = estimavarx(data, p, exo);

% Ajuste para betas_imp:
% betas tiene dimensiones (k * p + n_exo) filas.
% Debemos separar los coeficientes de las exógenas para las IRF.
n_exo_total = size(exo, 2);
betas_imp = betas(:, 1:end-n_exo_total); % Quita las ultimas columnas correspondientes a las exógenas

% Variance decomposition
for ii = 1:ne
    for rr = 1:ne
        eval(['resp' num2str(rr) 'imp' num2str(ii) ' = zeros(s, reps);']);
        eval(['var_resp' num2str(rr) 'imp' num2str(ii) ' = zeros(s, reps);']);
    end
end
%  Uhlig sample and IRFs building
i = 1;
while i <= reps
    % Draws for sigma
    R = randnm(0, inv(omega)/size(res,1), size(res,1));
    Sigma = inv(R * R');
    
    % Draws for B conditioned in Sigma
    xx = momentxx(data, p);
    B_cov = kron(inv(xx), Sigma);
    B = randnm(betas_imp(:), B_cov, 1);
    
    % Rearrange B
    ss = (1:k:length(B)+k)';
    storeb = zeros(k, length(ss)-1);
    for h = 1:length(ss)-1
        storeb(:,h) = B(ss(h):ss(h+1)-1,:);
    end
    
    % Cholesky and impulse vectors
        % Cholesky and impulse vectors
    C = transpose(chol(Sigma));   % Sigma = C * C'

    for ii = 1:ne
        v_imp = C(:, ii);         % choque de 1 s.d. estructural por defecto

        if USE_UNIT_SHOCK
            % Normalizar para que la respuesta contemporánea de la variable ii sea 1
            scale = v_imp(ii);
            if abs(scale) > 1e-12
                v_imp = v_imp / scale;
            else
                % Si algo sale raro (casi cero), mantenemos 1 s.d.
                warning('Shock %d: escala casi cero, se deja como 1 s.d.', ii);
            end
        end

        eval(['imp' num2str(ii) ' = v_imp;']);
    end

    
    % IRFs para cada shock
    for ii = 1:ne
        eval(['resp' num2str(ii) ' = impulse(storeb, imp' num2str(ii) ', s, p);' ]);
    end
    
    % Guardar IRFs por par
    for ii = 1:ne
        for rr = 1:ne
            eval(['resp' num2str(rr) 'imp' num2str(ii) '(:, i) = resp' num2str(ii) '(:, ' num2str(rr) ');' ]);
        end
    end
    
    % Variance decomposition
    for ii = 1:ne
        for rr = 1:ne
            aux = zeros(s,1);
            for rrr = 1:ne
                eval(['aux = aux + resp' num2str(rrr) '(:, ' num2str(rr) ').^2;']);
            end
            eval(['var_resp' num2str(rr) 'imp' num2str(ii) '(:, i) = (resp' num2str(ii) '(:, ' num2str(rr) ').^2) ./ aux;']);
        end
    end
    i = i + 1;
end
%% 6. IRFs
color1 = [0.5843, 0.8157, 0.9882];
color2 = [0, 0, 0.5];
mcol = 1;

% Adjust Time Axis for Monthly vs Quarterly
if monthly_data
    freq_divisor = 12; % 12 months per year
    xlabel_text = 'Years';
else
    freq_divisor = 4; % 4 quarters per year
    xlabel_text = 'Years';
end

x = (0:s-1)/freq_divisor; % X-axis values
years = 0:floor((s-1)/freq_divisor); % Ticks for integer years

variable_names = cellstr(var_order);
fprintf('\nVar variables in final order:\n');
for i = 1:numel(variable_names)
    fprintf('%d: %s\n', i, variable_names{i});
end
IMP = 0;
while isempty(IMP) || ~isnumeric(IMP) || IMP < 1 || IMP > k
    try
        IMP = input(sprintf('Choose the shock (1-%d): ', k));
        if isempty(IMP) || ~isnumeric(IMP) || IMP < 1 || IMP > k
            fprintf('Error: must be an int between 1 y %d.\n', k);
        end
    catch
        fprintf('Error: Invalid input.\n');
    end
end
fprintf('Plotting shock: %d: %s\n', IMP, variable_names{IMP});
figure('Name','IRFs VARX con Chebyshev exógeno','Color','w');
if use_precovid
    cutoff_text = sprintf('Corte en %s', CUTOFF);
else
    cutoff_text = 'Sin corte';
end
sgtitle(sprintf('IRF al shock de %s (Cholesky) | p=%d, m*=%d, s=%d, %s, %s, Exo=%s', ...
    variable_names{IMP}, p, m_star, s, cutoff_text,time_granularity,txt_exo), 'FontSize', 14);
for ii = 1:ne
    subplot(2, 4, ii);
    IR_data = eval(['resp' num2str(ii) 'imp' num2str(IMP)]);
    
    if use_95percent
        IR = quantile(IR_data', [0.025 0.50 0.975])'*100;
    else
        IR = quantile(IR_data', [0.16 0.50 0.84])'*100;
    end
    
    IRm = median(IR_data, 2) * 100;
    y1 = IR(:,1)';  
    y2 = IR(:,2)';  
    y3 = IR(:,3)';
    
    hold on;
    patch([x fliplr(x)], [y1 fliplr(y2)], color1, 'EdgeColor', 'none');
    patch([x fliplr(x)], [y2 fliplr(y3)], mcol*color1, 'EdgeColor', 'none');
    plot(x, zeros(1, s), 'r', 'LineWidth', 1);
    plot(x, IRm', 'Color', color2, 'LineWidth', 2);
    box on;
    title(variable_names{ii}, 'Interpreter', 'none');
    
    % Adjusted Axis Labels
    xticks(years);
    xticklabels(arrayfun(@num2str, years, 'UniformOutput', false));
    xlabel(xlabel_text);
    
    set(gca, 'Layer', 'top');
    hold off;
end
%%  7. Save final results
results_meta = struct();
results_meta.var_order = var_order;
results_meta.exo_names = external_exo_names;        % NUEVO: Nombres exógenas
results_meta.p = p;
results_meta.s = s;
results_meta.reps = reps;
results_meta.k = k;
results_meta.m = m;
results_meta.m_star = m_star;
results_meta.periods = periods_cma;
results_meta.shock_plotted_index = IMP;
results_meta.shock_plotted_name = variable_names{IMP};
save('varx_irf_results.mat', ...
     'betas','betas_imp','omega','res', ...
     'results_meta');
% (Opcional) Persistir IRFs y VDs por archivo
% for ii = 1:ne
%     for rr = 1:ne
%         eval(['save(sprintf(''irf_rr' num2str(rr) '_imp' num2str(ii) '.mat''), ''resp' num2str(rr) 'imp' num2str(ii) ''');']);
%         eval(['save(sprintf(''vd_rr'  num2str(rr) '_imp' num2str(ii) '.mat''), ''var_resp' num2str(rr) 'imp' num2str(ii) ''');']);
%     end
% end
toc;
fprintf('Saved as "varx_irf_results.mat".\n');