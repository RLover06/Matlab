function ajustePseudoVoigt_ZnO14()
    %% 1. Cargar espectro y corregir baseline
    % Leer los datos del espectro Raman (columna 1: desplazamiento Raman, columna 2: intensidad)
    data = readmatrix('C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Linea base\espectro_corregido.txt');    
    x = data(:,1);  % Eje X: desplazamiento Raman (cm^-1)
    y = data(:,2);  % Eje Y: intensidad del espectro

    % Corrección de la línea base utilizando algoritmo adaptativo
    bl = msbackadj(x, y, 'WindowSize', 50, 'StepSize', 25);
    y_corr = y - bl;  % Espectro corregido restando la línea base

    %% 2. Detección automática de picos
    % Buscar picos prominentes (mayores al 5% del máximo)
    [pks, locs] = findpeaks(y_corr, x, 'MinPeakProminence', 0.05*max(y_corr));
    N = numel(locs);  % Número de picos detectados

    % Nombres de modos vibracionales del ZnO conocidos
    modoNamesAll = { ...
        'E2(low)', 'TA', 'LA', 'TO', 'A1(TO)', 'E1(TO)', 'E2(high)', ...
        'LO', 'A1(LO)', 'E1(LO)', 'TO+TA', 'LO+LA', ...
        'E1(LO)+A1(LO)', 'TO+E2(high)', ...
        '2LA', '2TO', '2LO', 'E2(high)-E2(low)', ...
        'B1(high)-B1(low)', 'Surface Optical Phonon (SOP)'};

    % Centros teóricos de cada modo (cm^-1)
    modeCentersAll = [100, 150, 200, 400, 380, 410, 438, ...
                      580, 585, 590, 625, 700, 600, 550, ...
                      483, 980, 1105, 333, 284, 550];

    %% 3. Inicialización de parámetros para ajuste global (solo para visualizar)
    modos0 = [pks(:), locs(:), repmat([3, 3, 0.5], N, 1)];  % [Amplitud, Centro, Sigma, Gamma, Eta]
    b0 = [modos0(:); mean(y_corr)];  % Agregar constante de fondo al final

    % Definir límites para cada parámetro
    lb = []; ub = [];
    for i = 1:N
        A0 = modos0(i,1); x0 = modos0(i,2);
        lb = [lb; 0; x0-10; 0.5; 0.5; 0];       % Límites inferiores
        ub = [ub; 1.5*A0; x0+10; 20; 20; 1];    % Límites superiores
    end
    lb = [lb; -Inf]; ub = [ub; Inf];  % Constante de fondo sin límite

    %% 4. Ajuste global del espectro
    opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxFunctionEvaluations', 1e5);
    model = @(b,xdata) modeloPseudovoigt(b,xdata);
    bfit = lsqcurvefit(model, b0, x, y_corr, lb, ub, opts);

    %% 5. Gráfica del ajuste global
    yfit_corr = model(bfit, x);  % Ajuste sin baseline
    yfit = yfit_corr + bl;       % Se suma la baseline para comparar con datos originales

    figure('Units','normalized','Position',[0.1 0.1 0.8 0.7]); hold on;
    title('Ajuste de espectro Raman de ZnO con perfiles Pseudo-Voigt');
    xlabel('Desplazamiento Raman (cm^{-1})');
    ylabel('Intensidad (u.a.)');

    % Graficar espectro original y ajuste global
    h_exp = plot(x, y, 'k', 'LineWidth', 1.5, 'DisplayName', 'Experimental');
    h_fit = plot(x, yfit, 'r', 'LineWidth', 1.5, 'DisplayName', 'Fit Pseudo-Voigt');

    % Colores para cada marcador
    colors = lines(numel(modoNamesAll));

    % Agregar marcadores para cada modo vibracional
    handles = gobjects(numel(modoNamesAll), 1);
    for i = 1:numel(modoNamesAll)
        xc = modeCentersAll(i);
        yc = interp1(x, yfit, xc);
        handles(i) = plot(xc, yc, 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', modoNamesAll{i});
    text(xc, yc, modoNamesAll{i}, 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', 'Color', colors(i,:), ...
    'FontWeight', 'bold', 'Rotation', 60);

    end

    % Guardar gráfico como imagen
    exportgraphics(gcf, ...
      'C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Resultados\Ajuste_PseudoVoigt_ZnO.png', ...
      'Resolution', 300);

    %% 6. Ajuste local por cada modo vibracional
    resultados = cell(numel(modoNamesAll), 7);  % Inicializar tabla de resultados
    for i = 1:numel(modoNamesAll)
        nombre = modoNamesAll{i};
        x0_teorico = modeCentersAll(i);
        resultados{i,1} = nombre;
        resultados{i,2} = x0_teorico;

        % Definir ventana de ajuste alrededor del centro teórico
        rango = 25;
        mascara = x >= (x0_teorico - rango) & x <= (x0_teorico + rango);
        x_sub = x(mascara);
        y_sub = y_corr(mascara);

        if numel(x_sub) < 10
            resultados{i,3:7} = {NaN, NaN, NaN, NaN, NaN};
            continue;
        end

        % Parámetros iniciales del ajuste local
        A0 = max(y_sub);
        sigma0 = 5; gamma0 = 5; eta0 = 0.5;
        p0 = [A0, x0_teorico, sigma0, gamma0, eta0];
        lb = [0, x0_teorico - 10, 0.5, 0.5, 0];
        ub = [2*A0, x0_teorico + 10, 20, 20, 1];

        % Ajuste con lsqcurvefit al pico individual
        model_pv = @(p, xdata) p(1) * pseudoVoigt(xdata, p(2), p(3), p(4), p(5));
        opts = optimoptions('lsqcurvefit','Display','off');
        try
            pfit = lsqcurvefit(model_pv, p0, x_sub, y_sub, lb, ub, opts);
            resultados{i,3} = pfit(2); % Centro ajustado
            resultados{i,4} = pfit(3); % Sigma
            resultados{i,5} = pfit(4); % Gamma
            resultados{i,6} = pfit(5); % Eta
            resultados{i,7} = pfit(1); % Amplitud
        catch
            resultados{i,3:7} = {NaN, NaN, NaN, NaN, NaN};
        end
    end

    % Exportar resultados a archivo Excel
    T = cell2table(resultados, 'VariableNames', {
        'ModoNormal', 'CentroTeorico_cm1', 'CentroAjustado_cm1', ...
        'Sigma', 'Gamma', 'Eta', 'Amplitud'});
    writetable(T, 'ajuste_pseudovoigt_local.xlsx');
    disp('Archivo Excel con parámetros ajustados exportado correctamente.');
end

%% Modelo Pseudo-Voigt para ajuste global
function y = modeloPseudovoigt(b,x)
    N = (length(b)-1)/5;
    y = zeros(size(x));
    for j = 1:N
        idx = (j-1)*5 + (1:5);
        A = b(idx(1)); x0 = b(idx(2));
        s = b(idx(3)); g = b(idx(4)); eta = b(idx(5));
        y = y + A * pseudoVoigt(x, x0, s, g, eta);
    end
    y = y + b(end);  % Agregar constante de fondo
end

%% Función del perfil Pseudo-Voigt normalizado
function V = pseudoVoigt(x, x0, sigma, gamma, eta)
    G = exp(-((x - x0).^2)/(2*sigma^2)) / (sigma * sqrt(2*pi));  % Gaussiana
    L = (gamma/pi) ./ ((x - x0).^2 + gamma^2);                   % Lorentziana
    V = eta * L + (1 - eta) * G;                                 % Mezcla ponderada
    V = V / max(V);  % Normalizar el pico a 1
end
