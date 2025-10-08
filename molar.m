% Datos
x = [0, 0.01, 0.02, 0.03, 0.04, 0.05]; % Eje X
II = [2.297423524, 0.959289016, 1.536413285, 1.072074826, 1.050525176, 0.406505807]; % Eje Y

% Crear la figura
figure;
plot(x, II, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.1, 0.4, 0.7], 'MarkerEdgeColor', 'k');
grid on;

% Etiquetas
xlabel('Relación molar del dopado x', 'FontSize', 12);
ylabel('Razón de intensidades I(LO)/ I(E_2^{high}) ', 'FontSize', 12);
title('Relación molar x  vs  I(LO)/ I(E_2^{high})', 'FontSize', 14);

% Mejoras visuales
xlim([0 0.06]);
ylim([0 2.5]);

