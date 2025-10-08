function ajustePseudoVoigt_ZCe3()
    %% 1. Cargar espectro y corregir baseline
    % Leer datos desde archivo (columna 1 = eje X: desplazamiento Raman, columna 2 = eje Y: intensidad)
    data = readmatrix('C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Linea base\espectro_corregido.txt');    
    x = data(:,1);  % Extraer eje X (desplazamiento Raman en cm^-1)
    y = data(:,2);  % Extraer eje Y (intensidad)

    % Aplicar corrección de línea base usando método adaptativo
    bl = msbackadj(x, y, 'WindowSize', 50, 'StepSize', 25); 
    y_corr = y - bl;  % Espectro corregido restando la línea base

    %% 2. Detección automática de picos
    % Buscar picos con prominencia mayor al 5% del valor máximo
    [pks, locs] = findpeaks(y_corr, x, 'MinPeakProminence', 0.05*max(y_corr));
    N = numel(locs);  % Número total de picos detectados

    % Lista de nombres de modos vibracionales teóricos del ZnO
    modoNamesAll = { ...
        'E2(low)', 'TA', 'LA', 'TO', 'A1(TO)', 'E1(TO)', 'E2(high)', ...
        'LO', 'A1(LO)', 'E1(LO)', 'TO+TA', ...
        'TO+E2(high)', ...
        '2LA; 2B1(low)', '2TO', '2LO', 'E2(high)-E2(low)', ...
        'B1(high)-B1(low)', '2TA; 2E2(low)', 'LA + LO','Ce0_2'};

    % Centros teóricos (en cm^-1) de cada modo vibracional
    modeCentersAll = [102, 131, 270, 404, 380, 412.2, 440, ...
                      562.8, 574.9, 590, 620, 350, ...
                      537.7, 983.6, 1111.5, 334, 284, 207, 813.8, 466];

    % FWHM teóricos para comparación (en cm^-1), algunos son NaN si no se conocen
    FWHM_teorico = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]';

    %% 3. Inicialización de parámetros para ajuste global
    % Crear arreglo inicial: [Amplitud, Centro, Sigma, Gamma, Eta]
    modos0 = [pks(:), locs(:), repmat([3, 3, 0.5], N, 1)];
    b0 = [modos0(:); mean(y_corr)];  % Vector de parámetros iniciales + fondo constante

    % Establecer límites inferiores y superiores para cada parámetro
    lb = []; ub = [];
    for i = 1:N
        A0 = modos0(i,1); x0 = modos0(i,2);
        lb = [lb; 0; x0-10; 0.5; 0.5; 0];       % Límite inferior
        ub = [ub; 1.5*A0; x0+10; 20; 20; 1];    % Límite superior
    end
    lb = [lb; -Inf]; ub = [ub; Inf];  % Sin límite para la constante de fondo

    %% 4. Ajuste global del espectro
    opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxFunctionEvaluations', 1e5);  % Configuración del optimizador
    model = @(b,xdata) modeloPseudovoigt(b,xdata);  % Definir función objetivo
    bfit = lsqcurvefit(model, b0, x, y_corr, lb, ub, opts);  % Ejecutar el ajuste global

    %% 5. Graficar ajuste global con resultados
    yfit_corr = model(bfit, x);     % Curva ajustada sin línea base
    yfit = yfit_corr + bl;          % Se suma la línea base para comparación con el espectro original

    figure('Units','normalized','Position',[0.1 0.1 0.8 0.7]); hold on;
    title('Ajuste de espectro Raman de ZCe3 con perfiles Pseudo-Voigt');
    xlabel('Desplazamiento Raman (cm^{-1})');
    ylabel('Intensidad (u.a.)');

    % Graficar espectro experimental y ajuste
    h_exp = plot(x, y, 'k', 'LineWidth', 1.5, 'DisplayName', 'Experimental');
    h_fit = plot(x, yfit, 'r', 'LineWidth', 1.5, 'DisplayName', 'Fit Pseudo-Voigt');

    % Generar colores y marcadores para cada modo
    colors = lines(numel(modoNamesAll));
    handles = gobjects(numel(modoNamesAll), 1);
    for i = 1:numel(modoNamesAll)
        xc = modeCentersAll(i);                      % Centro teórico
        yc = interp1(x, yfit, xc);                   % Interpolar intensidad en ese punto
        handles(i) = plot(xc, yc, 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', modoNamesAll{i});         % Añadir marcador
        text(xc, yc, modoNamesAll{i}, 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', 'Color', colors(i,:), ...
            'FontWeight', 'bold');                   % Etiqueta del modo vibracional
    end

    legend([h_exp, h_fit, handles(:)'], 'Location', 'northeastoutside');  % Leyenda
    grid on;

    % Guardar figura como archivo PNG
    exportgraphics(gcf, ...
      'C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Resultados\Ajuste_PseudoVoigt_ZnO.png', ...
      'Resolution', 300);

    %% 6. Ajuste local por cada modo vibracional y cálculo de parámetros físicos
    resultados = cell(numel(modoNamesAll), 9);  % Inicializar tabla de resultados
    for i = 1:numel(modoNamesAll)
        nombre = modoNamesAll{i};
        x0_teorico = modeCentersAll(i);
        resultados{i,1} = nombre;
        resultados{i,2} = x0_teorico;

        % Extraer una ventana alrededor del modo teórico
        rango = 25;
        mascara = x >= (x0_teorico - rango) & x <= (x0_teorico + rango);
        x_sub = x(mascara);
        y_sub = y_corr(mascara);

        if numel(x_sub) < 10
            resultados{i,3:9} = {NaN, NaN, NaN, NaN, NaN, NaN, NaN};
            continue;
        end

        % Estimar sigma y gamma a partir del FWHM teórico si está disponible
        if ~isnan(FWHM_teorico(i))
            sigma0 = FWHM_teorico(i)/2.3548;
            gamma0 = FWHM_teorico(i)/2;
            lb = [0, x0_teorico - 10, sigma0*0.5, gamma0*0.5, 0];
            ub = [2*max(y_sub), x0_teorico + 10, sigma0*2, gamma0*2, 1];
        else
            sigma0 = 5; gamma0 = 5;
            lb = [0, x0_teorico - 10, 0.5, 0.5, 0];
            ub = [2*max(y_sub), x0_teorico + 10, 20, 20, 1];
        end

        A0 = max(y_sub); eta0 = 0.5;
        p0 = [A0, x0_teorico, sigma0, gamma0, eta0];

        % Ajustar pico individual con perfil Pseudo-Voigt
        model_pv = @(p, xdata) p(1) * pseudoVoigt(xdata, p(2), p(3), p(4), p(5));
        opts = optimoptions('lsqcurvefit','Display','off');
        try
            pfit = lsqcurvefit(model_pv, p0, x_sub, y_sub, lb, ub, opts);
            A = pfit(1); xc = pfit(2); sig = pfit(3); gam = pfit(4); eta = pfit(5);

            % Cálculo del área bajo la curva ajustada
            y_model = A * pseudoVoigt(x_sub, xc, sig, gam, eta);
            area = trapz(x_sub, y_model);

            % Cálculo de FWHM ajustado
            FWHM_adj = eta * 2 * gam + (1 - eta) * 2.3548 * sig;

            % Almacenar resultados
            resultados(i,3:9) = {xc, sig, gam, eta, A, area, FWHM_adj};
        catch
            resultados{i,3:9} = {NaN, NaN, NaN, NaN, NaN, NaN, NaN};
        end
    end

    %% 7. Exportar tabla final con todos los parámetros y comparación
    T = cell2table(resultados, 'VariableNames', {
        'ModoNormal', 'CentroTeorico_cm1', 'CentroAjustado_cm1', ...
        'Sigma', 'Gamma', 'Eta', 'Amplitud', 'AreaBajoCurva', 'FWHM_Ajustado_cm1'});
    T.FWHM_Teorico_cm1 = FWHM_teorico;  % Agregar columna de FWHM teórico
    writetable(T, 'ajuste_pseudovoigt_local.xlsx');  % Guardar como Excel

    disp('Archivo Excel con parámetros ajustados y comparaciones exportado correctamente.');
end

%% Modelo para ajuste global (suma de picos Pseudo-Voigt + fondo)
function y = modeloPseudovoigt(b,x)
    N = (length(b)-1)/5;
    y = zeros(size(x));
    for j = 1:N
        idx = (j-1)*5 + (1:5);
        A = b(idx(1)); x0 = b(idx(2));
        s = b(idx(3)); g = b(idx(4)); eta = b(idx(5));
        y = y + A * pseudoVoigt(x, x0, s, g, eta);  % Sumar cada pico
    end
    y = y + b(end);  % Agregar la constante de fondo
end

%% Perfil Pseudo-Voigt normalizado
function V = pseudoVoigt(x, x0, sigma, gamma, eta)
    G = exp(-((x - x0).^2)/(2*sigma^2)) / (sigma * sqrt(2*pi));  % Parte gaussiana
    L = (gamma/pi) ./ ((x - x0).^2 + gamma^2);                   % Parte lorentziana
    V = eta * L + (1 - eta) * G;                                 % Combinación L-G
    V = V / max(V);  % Normalizar a valor máximo de 1
end
