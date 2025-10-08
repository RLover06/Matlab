function ajustePseudoVoigt_ZC1()
    %% 1. Cargar espectro y corregir baseline
    data = readmatrix('C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Linea base\espectro_corregido.txt');    
    x = data(:,1); y = data(:,2);
    bl = msbackadj(x, y, 'WindowSize', 50, 'StepSize', 25);
    y_corr = y - bl;

    %% 2. Detección automática de picos
    [pks, locs] = findpeaks(y_corr, x, 'MinPeakProminence', 0.05*max(y_corr));
    N = numel(locs);

    modoNamesAll = { ...
        'E2(low)', 'TA', 'LA', 'TO', 'A1(TO)', 'E1(TO)', 'E2(high)', ...
        'LO', 'A1(LO)', 'E1(LO)', 'TO+TA', ...
        'TO+E2(high)', ...
        '2B1(low)', '2TO', '2LO', 'E2(high)-E2(low)', ...
        'B1(high)-B1(low)', '2E2(low)', 'LA + LO'};

    modeCentersAll = [102, 131, 270, 404, 380, 412.2, 440, ...
                      562.8, 574.9, 590, 620, 350, ...
                      537.7, 983.6, 1111.5, 334, 284, 207, 813.8];

      % FWHM teóricos para comparación (en cm^-1), algunos son NaN si no se conocen
    FWHM_teorico = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]';

    %% 3. Inicialización para ajuste global
    modos0 = [pks(:), locs(:), repmat([3, 3, 0.5], N, 1)];
    b0 = [modos0(:); mean(y_corr)];

    lb = []; ub = [];
    for i = 1:N
        A0 = modos0(i,1); x0 = modos0(i,2);
        lb = [lb; 0; x0-10; 0.5; 0.5; 0];
        ub = [ub; 1.5*A0; x0+10; 20; 20; 1];
    end
    lb = [lb; -Inf]; ub = [ub; Inf];

    %% 4. Ajuste global
    opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxFunctionEvaluations', 1e5);
    model = @(b,xdata) modeloPseudovoigt(b,xdata);
    bfit = lsqcurvefit(model, b0, x, y_corr, lb, ub, opts);

    %% 5. Gráfica
    yfit_corr = model(bfit, x);
    yfit = yfit_corr + bl;

    figure('Units','normalized','Position',[0.1 0.1 0.8 0.7]); hold on;
   % title('Ajuste de espectro Raman de ZCe5 con perfiles Pseudo-Voigt');
    xlabel('Desplazamiento Raman (cm^{-1})'); ylabel('Intensidad (u.a.)');

    h_exp = plot(x, y, 'k', 'LineWidth', 1.5, 'DisplayName', 'ZnO');
    h_fit = plot(x, yfit, 'r', 'LineWidth', 1.5, 'DisplayName', 'Fit Pseudo-Voigt');

    colors = lines(numel(modoNamesAll));
    handles = gobjects(numel(modoNamesAll), 1);
    for i = 1:numel(modoNamesAll)
        xc = modeCentersAll(i);
        yc = interp1(x, yfit, xc);
        handles(i) = plot(xc, yc, 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', modoNamesAll{i});
        text(xc, yc, modoNamesAll{i}, 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', 'Color', colors(i,:), ...
            'FontWeight', 'bold');
    end
    legend([h_exp, h_fit, handles(:)'], 'Location', 'northeastoutside');
    grid on;

    exportgraphics(gcf, ...
      'C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Resultados\Ajuste_PseudoVoigt_ZnO.png', ...
      'Resolution', 300);

    %% 6. Ajuste local con restricciones
    resultados = cell(numel(modoNamesAll), 10);
    for i = 1:numel(modoNamesAll)
        nombre = modoNamesAll{i};
        x0_teorico = modeCentersAll(i);
        resultados{i,1} = nombre;
        resultados{i,2} = x0_teorico;

        rango = 25;
        mascara = x >= (x0_teorico - rango) & x <= (x0_teorico + rango);
        x_sub = x(mascara);
        y_sub = y_corr(mascara);

        if numel(x_sub) < 10
            resultados(i,3:10) = {NaN};
            continue;
        end

        A0 = max(y_sub);
        sigma0 = 5; gamma0 = 5; eta0 = 0.5;
        p0 = [A0, x0_teorico, sigma0, gamma0, eta0];

        if ~isnan(FWHM_teorico(i))
            fwhm_target = FWHM_teorico(i);
            sigma_min = max(0.5, 0.3 * fwhm_target);
            sigma_max = 2 * fwhm_target;
            gamma_min = max(0.5, 0.3 * fwhm_target);
            gamma_max = 2 * fwhm_target;
        else
            sigma_min = 0.5; sigma_max = 20;
            gamma_min = 0.5; gamma_max = 20;
        end

        lb = [0, x0_teorico - 10, sigma_min, gamma_min, 0];
        ub = [2*A0, x0_teorico + 10, sigma_max, gamma_max, 1];

        model_pv = @(p, xdata) p(1) * pseudoVoigt(xdata, p(2), p(3), p(4), p(5));
        try
            pfit = lsqcurvefit(model_pv, p0, x_sub, y_sub, lb, ub, opts);
            A = pfit(1); x0 = pfit(2); sigma = pfit(3); gamma = pfit(4); eta = pfit(5);
            FWHM = (1 - eta) * 2.355 * sigma + eta * 2 * gamma;
            dx = mean(diff(x_sub));
            y_model = model_pv(pfit, x_sub);
            area = trapz(x_sub, y_model);

            resultados(i,3:10) = {x0, sigma, gamma, eta, A, area, FWHM, FWHM_teorico(i)};
        catch
            resultados(i,3:10) = {NaN};
        end
    end

    %% 7. Exportar a Excel
    T = cell2table(resultados, 'VariableNames', {
        'ModoNormal', 'CentroTeorico_cm1', 'CentroAjustado_cm1', ...
        'Sigma', 'Gamma', 'Eta', 'Amplitud', ...
        'AreaBajoCurva', 'FWHM_Ajustado_cm1', 'FWHM_Teorico_cm1'});

    writetable(T, 'ajuste_pseudovoigt_local.xlsx');
    disp('Archivo Excel con parámetros ajustados exportado correctamente.');
end

%% Modelo Pseudo-Voigt global
function y = modeloPseudovoigt(b,x)
    N = (length(b)-1)/5;
    y = zeros(size(x));
    for j = 1:N
        idx = (j-1)*5 + (1:5);
        A = b(idx(1)); x0 = b(idx(2));
        s = b(idx(3)); g = b(idx(4)); eta = b(idx(5));
        y = y + A * pseudoVoigt(x, x0, s, g, eta);
    end
    y = y + b(end);
end

%% Perfil Pseudo-Voigt normalizado
function V = pseudoVoigt(x, x0, sigma, gamma, eta)
    G = exp(-((x - x0).^2)/(2*sigma^2)) / (sigma * sqrt(2*pi));
    L = (gamma/pi) ./ ((x - x0).^2 + gamma^2);
    V = eta * L + (1 - eta) * G;
    V = V / max(V);
end
