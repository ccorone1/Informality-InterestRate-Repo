%% Data processing and SVAR with Chevysheb estimation with input selection
%
% By Carlos Coronel (Last update: 06 february, 2026)                                                  
%
% CORRECCIÓN ACTUALIZADA:
% 1. Calcula CMA sobre Full Sample (evita end-point bias).
% 2. Corta la muestra.
% 3. Calcula Chebyshev sobre la muestra cortada (garantiza dominio [-1, 1]).
%

clear all; clc; close all; tic;
%addpath("C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\code\matlab");
%addpath("C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\code\matlab\20_11_25");
addpath("C:\Users\Carlos Coronel\Documents\Informalidad\Informality-InterestRate\Code\calls");


%%  1. Parameters
% default index
default_var_indices = [1,2,3,4];
USE_INTERACTIVE_SELECTION = true; 
use_precovid = false; %crop to cutoff
use_cosine_def = true;
use_68percent = true; %You can choose to plot one or both of the CI, if both are false none will be plotted
use_95percent = true; 
monthly_data = true; %Usar datos mensuales
PLOT_CHEB_TREND_SUBPLOTS = true;   % true/false: figura adicional con series + tendencia Chebyshev 
USE_UNIT_SHOCK = true;

% --- Parámetros del Modelo ---
if monthly_data
    CUTOFF = "12/01/2019"; 
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\mensuales_2005.xlsx";
    time_granularity = "Mensual";
else
    CUTOFF = "2019Q4";
    filename = "C:\Users\Carlos Coronel\Documents\Respaldo Memoria\leyva_research\data\clean\trimestrales_primero.xlsx";
    time_granularity = "Trimestral";
end
m = 10;      % Chebyshev max order
m_star = 6; % Degree of the polynomial to be used
p = 1;      % VAR lags
s = 49;     % IRF horizon (steps)
reps = 1000;  % Monte carlo repetitions

%% 2. Loading, selection and processing
%  loading
data = readtable(filename);
% Formatting periods
if monthly_data
    try
        if isdatetime(data{:,1})
             periods = string(datetime(data{:,1}, 'Format', 'MM/dd/yyyy'));
        else
             periods = string(data{:,1});
        end
    catch
        warning('Hubo un problema formateando las fechas. Se usarán tal cual vienen.');
        periods = string(data{:,1});
    end
else
    periods = string(data{:,1});    
end

Ytbl_raw = data(:,2:end);       
allVarNames = Ytbl_raw.Properties.VariableNames;

% --- Dynamic selection (Endógenas) ---
if USE_INTERACTIVE_SELECTION 
    fprintf('Data Columns:\n');
    for i = 1:numel(allVarNames)
        fprintf('%d: %s\n', i, allVarNames{i});
    end
    sel_idx = input('Indices para ENDÓGENAS (VAR / Cholesky Order) (ej: [1, 3]): ');
    
    % --- Selección de Exógenas ---
    fprintf('\nSelecciona las variables EXÓGENAS adicionales (además del polinomio).\n');
    fprintf('Deja vacío y pulsa Enter si no quieres ninguna.\n');
    sel_exo_idx = input('Indices para EXÓGENAS (ej: [5, 6] o []): ');
else
    sel_idx = default_var_indices;
    sel_exo_idx = []; 
end

% Filtrar Endógenas
Ytbl = Ytbl_raw(:, sel_idx);
var_order = Ytbl.Properties.VariableNames;

% Filtrar Exógenas (Raw)
if ~isempty(sel_exo_idx)
    Ytbl_exo_raw = Ytbl_raw(:, sel_exo_idx);
else
    Ytbl_exo_raw = table(); 
end
fprintf('Var variables in final order:\n');
disp(var_order');

% NOTA: El corte de muestra se ha movido a la Sección 3 para permitir que el CMA use datos futuros.

[n, p_selected] = size(Ytbl); 

% Extraer nombres de las exógenas externas ---
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

%% 3. CMA & CHEBYSHEV (CORREGIDO: CÁLCULO -> CORTE -> POLINOMIO)

% 3.1 CÁLCULO DE CMA (Sobre muestra completa para evitar sesgo de borde)
if monthly_data
    if n < 13
        error('Error: Not enough observations for Monthly CMA (Need > 12).');
    end
    idx_cma = 7:(n-6);
    cma_func = @cma_m; 
else
    if n < 5
        error('Error: Not enough observations for Quarterly CMA (Need > 4).');
    end
    idx_cma = 3:(n-2);
    cma_func = @cma_q; 
end
n_cma = numel(idx_cma);
cma_series = NaN(n_cma, p_selected);

% Execute CMA for endo
for j = 1:p_selected
    y = double(Ytbl{:,j});
    try
        y_cma = cma_func(y); 
        if length(y_cma) ~= n_cma
            error('CMA output length mismatch.');
        end
        cma_series(:,j) = y_cma;
    catch ME
        warning('Error calculating CMA for column %s: %s', var_order{j}, ME.message);
    end
end

% Execute CMA for exo
cma_exo_series = [];
if ~isempty(sel_exo_idx)
    n_exo_selected = length(sel_exo_idx);
    cma_exo_series = NaN(n_cma, n_exo_selected);
    for j = 1:n_exo_selected
        y = double(Ytbl_exo_raw{:,j}); 
        try
            y_cma = cma_func(y); 
            cma_exo_series(:,j) = y_cma;
        catch ME
            warning('Error en CMA exógena %d: %s', j, ME.message);
        end
    end
end

endo_raw = cma_series;
periods_cma = periods(idx_cma);

% 3.2 LÓGICA DE CORTE 
% Cortamos los datos YA SUAVIZADOS pero ANTES de calcular el polinomio
if use_precovid
    fprintf('\n--- Buscando fecha de corte: %s ---\n', CUTOFF);
    cutoff_idx = find(strcmp(periods_cma, CUTOFF));
    
    if isempty(cutoff_idx)
        warning('La fecha "%s" no se encontró. Se usará muestra completa.', CUTOFF);
    else
        % Cortamos las series endógenas y exógenas (externas)
        endo_raw = endo_raw(1:cutoff_idx, :);
        periods_cma = periods_cma(1:cutoff_idx);
        if ~isempty(cma_exo_series)
            cma_exo_series = cma_exo_series(1:cutoff_idx, :);
        end
        
        fprintf('Datos cortados exitosamente en %s. Nueva longitud: %d obs.\n', CUTOFF, cutoff_idx);
        % Actualizamos n_cma para que el polinomio se genere con el largo correcto
        n_cma = cutoff_idx; 
    end
end

%% 3.3 GENERACIÓN DEL POLINOMIO (nodos Chebyshev / DCT)
% Nota: Esta base es por índice de tiempo (1..n_cma), no por calendario.
% Requiere que existan en el path: chevyp.m, chevyp_ort0.m, chevyp_ort1.m, chevyp_ort2.m, chevyp_ort3.m

tt = (1:n_cma)';   % índice temporal

su = [];
Jmax = ceil(m_star/2);

for j = 1:Jmax
    % Término impar (orden 2j-1) ortogonalizado a [1, t/n] (y a impares previos)
    if (2*j - 1) <= m_star
        su = [su, chevyp_ort2(j, n_cma, tt)];
    end

    % Término par (orden 2j)
    if (2*j) <= m_star
        su = [su, chevyp_ort3(j, n_cma, tt)];
    end
end

% Incluye el término lineal "ortogonal" (sin constante)
T_cheb = [chevyp_ort1(n_cma, tt), su];

% (Recomendado) Centrar para asegurar ortogonalidad con el intercept del VARX
T_cheb = T_cheb - mean(T_cheb, 1);

% Conteo real de términos que entran como exógenos (incluye el lineal)
n_cheb_terms = size(T_cheb, 2);
fprintf('Chebyshev terms in exo (incl. linear): %d | max order k<=%d\n', n_cheb_terms, m_star);


%% 4. log(Variable) y Construcción Final
% Log a las endógenas 
endo_log = log(endo_raw); 

% Log a las exógenas externas (si existen)
exo_external_final = [];
if ~isempty(cma_exo_series)
    if any(cma_exo_series(:) <= 0)
        warning('Cuidado: Hay valores <= 0 en las exógenas externas antes del log.');
    end
    exo_external_final = log(cma_exo_series);
end

% Input variables final
endo = double(endo_log);
% CONCATENACIÓN: [Polinomio, Series_Externas]
exo = [double(T_cheb), double(exo_external_final)];

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

fprintf('Estimando VARX con %d variables exógenas (%d Chebyshev incl. linear + %d Externas)\n', ...
    size(exo,2), size(T_cheb,2), size(exo_external_final,2));
% VARX con OLS 
[betas, xxx, omega, res] = estimavarx(data, p, exo);



% Ajuste para betas_imp:
n_exo_total = size(exo, 2);
betas_imp = betas(:, 1:end-n_exo_total);

%%%%%%%%%

% ---- Correct X'X for [const+lags] conditional on exo using Schur complement ----
k = size(data,2);
r = 1 + k*p;               % number of coefficients per equation excluding exo
q = size(xxx,1);           % total regressors = r + n_exo_total

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

    xx_cond = GAM11 - GAM12 * GAM22inv * GAM12';   % = Z1 * Mx * Z1'
else
    xx_cond = xxx; % no exo case
end

%%%%%%%%%
 
assert(size(betas,2) == 1 + k*p + n_exo_total, 'betas columns do not match [const + lags + exo].');
assert(all(size(xx_cond) == [1+k*p, 1+k*p]), 'xx_cond must be r-by-r.');

% Variance decomposition init
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
    %xx = momentxx(data, p); replaced with xxx above
    B_cov = kron(inv(xx_cond), Sigma);
    B = randnm(betas_imp(:), B_cov, 1);
    
    % Rearrange B
    ss = (1:k:length(B)+k)';
    storeb = zeros(k, length(ss)-1);
    for h = 1:length(ss)-1
        storeb(:,h) = B(ss(h):ss(h+1)-1,:);
    end
    
    % Cholesky and impulse vectors
    C = transpose(chol(Sigma));   % Sigma = C * C'
    for ii = 1:ne
        v_imp = C(:, ii);         % choque de 1 s.d. estructural por defecto
        if USE_UNIT_SHOCK
            scale = v_imp(ii);
            if abs(scale) > 1e-12
                v_imp = v_imp / scale;
            else
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

%% 6. IRFs (INTEGRADO: selección 68/95 con true/false y colores distintos)

% --- Colores ---
%col95 = [0.5843, 0.8157, 0.9882];   % 95% azul claro
col95 = [0.25, 0.55, 0.85];   % más oscuro que el original

col68 = [0.0000, 0.0000, 0.5000];   % 68% azul oscuro
color_med = [0, 0, 0.5];            % línea mediana (puedes poner [0 0 0] si prefieres)
mcol = 1;

% --- Eje X ---
if monthly_data
    freq_divisor = 12;
    xlabel_text  = 'Years';
else
    freq_divisor = 4;
    xlabel_text  = 'Years';
end
x = (0:s-1)/freq_divisor;
years = 0:floor((s-1)/freq_divisor);

% --- Nombres variables ---
variable_names = cellstr(var_order);

fprintf('\nVar variables in final order:\n');
for i = 1:numel(variable_names)
    fprintf('%d: %s\n', i, variable_names{i});
end

% --- Elegir shock a graficar ---
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

% --- Figura / título ---
figure('Name','IRFs VARX con Chebyshev exógeno','Color','w');

if use_precovid
    cutoff_text = sprintf('Sin covid, corte en %s', CUTOFF);
else
    cutoff_text = 'Con covid';
end

sgtitle(sprintf('IRF al shock de %s (Cholesky) | p=%d, ChebTerms=%d (k<=%d), s=%d, %s, %s, Exo=%s', ...
    variable_names{IMP}, p, size(T_cheb,2), m_star, s, cutoff_text, time_granularity, txt_exo), 'FontSize', 14);

% --- Layout de subplots (robusto si ne cambia) ---
ncols = 4;
nrows = 2;

for ii = 1:(nrows*ncols)
    subplot(nrows, ncols, ii);

    % Si no hay variable para este slot, lo apagas y sigues
    if ii > ne
        axis off;
        continue
    end

    IR_data = eval(['resp' num2str(ii) 'imp' num2str(IMP)]);   % (s x reps)
    IRm     = median(IR_data, 2);                               % (s x 1)

    hold on;

    % 95% primero (fondo), luego 68% encima
    if use_95percent
        IR_95 = quantile(IR_data', [0.025 0.50 0.975])';  % (s x 3)
        lo95 = IR_95(:,1)';  hi95 = IR_95(:,3)';

        patch([x fliplr(x)], [lo95 fliplr(hi95)], col95, ...
              'EdgeColor','none', 'FaceAlpha', 0.25);
    end

    if use_68percent
        IR_68 = quantile(IR_data', [0.16 0.50 0.84])';     % (s x 3)
        lo68 = IR_68(:,1)';  hi68 = IR_68(:,3)';

        patch([x fliplr(x)], [lo68 fliplr(hi68)], col68, ...
              'EdgeColor','none', 'FaceAlpha', 0.35);
    end

    % Línea cero y mediana
    plot(x, zeros(1, s), 'r', 'LineWidth', 1);
    plot(x, IRm', 'Color', color_med, 'LineWidth', 2);

    box on;
    title(variable_names{ii}, 'Interpreter', 'none');

    xticks(years);
    xticklabels(arrayfun(@num2str, years, 'UniformOutput', false));
    xlabel(xlabel_text);

    set(gca, 'Layer', 'top');
    hold off;
end


%% 6.1 Subplots: Series vs tendencia Chebyshev
if PLOT_CHEB_TREND_SUBPLOTS

    if monthly_data
        try
            xdt = datetime(periods_cma, "InputFormat","MM/dd/yyyy");
        catch
            xdt = datetime(periods_cma);
        end
    else
        try
            xdt = datetime(periods_cma, "InputFormat","uuuu'Q'Q");
        catch
            xdt = datetime(periods_cma);
        end
        xdt = dateshift(xdt, 'end', 'quarter');
    end

    % --- Datos: endo ya está en log(CMA) y cortado ---
    Yall = endo;                % (T x ne)
    Tobs = size(Yall,1);
    tt   = (1:Tobs)';

    % --- Construye la base EXACTA para tendencia (incluye constante) ---
    % X = [chevyp_ort0, chevyp_ort1, su], con su formado por (impar ort2, par ort3)
    su = [];
    Jmax = ceil(m_star/2);
    for j = 1:Jmax
        if (2*j - 1) <= m_star
            su = [su, chevyp_ort2(j, Tobs, tt)];
        end
        if (2*j) <= m_star
            su = [su, chevyp_ort3(j, Tobs, tt)];
        end
    end

    Xtrend = [chevyp_ort0(Tobs), chevyp_ort1(Tobs, tt), su];   % (T x Ktrend)

    % --- Figura aparte ---
    figure('Name','Series vs Tendencia Chebyshev','Color','w');

    ncols_tr = 2;                       % más legible que 4 para series
    nrows_tr = ceil(ne / ncols_tr);

    for j = 1:(nrows_tr*ncols_tr)
        subplot(nrows_tr, ncols_tr, j);

        if j > ne
            axis off;
            continue
        end

        Y = Yall(:, j);

        A = (Xtrend' * Xtrend);
        b = (Xtrend' * Y);
        beta_ols = gaussj(A, b);

        YS = Xtrend * beta_ols;   % tendencia Cheb

        plot(xdt, Y,  'LineWidth', 1.0); hold on;
        plot(xdt, YS, 'LineWidth', 2.0);
        grid on; box on;

        title(variable_names{j}, 'Interpreter','none');

        % Etiquetas: solo abajo para no saturar
        if j > (nrows_tr-1)*ncols_tr
            xlabel('Time');
        end

        % Leyenda solo en el primer subplot
        if j == 1
            legend('Y (log, CMA)','Tendencia Chebyshev','Location','best');
        end

        hold off;
    end

    sgtitle(sprintf('Series vs Tendencia Chebyshev | m*=%d | %s | %s', ...
        m_star, cutoff_text, time_granularity), 'FontSize', 14);

end



%%  7. Save final results
results_meta = struct();
results_meta.var_order = var_order;
results_meta.exo_names = external_exo_names;        
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

toc;
fprintf('Saved as "varx_irf_results.mat".\n');