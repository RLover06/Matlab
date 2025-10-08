% Establecer la carpeta donde se encuentra el archivo
carpeta = 'C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Resultados';
cd(carpeta); % Cambia el directorio de trabajo actual a esa carpeta

% Verificar si el archivo existe
if ~exist('ZCe3.txt', 'file') % Comprueba si el archivo existe en la carpeta
    error('El archivo ZCe3.txt no se encuentra en la carpeta especificada.');
    % Muestra error si no existe
end
%==============================================================================
% Cargar datos del espectro Raman
data = load('ZCe3.txt'); % Carga los datos del archivo como matriz numérica
raman_shift = data(:,1);  % Columna 1: eje X (desplazamiento Raman)
intensity = data(:,2);    % Columna 2: eje Y (intensidad)
%===============================================================================
% --- Escalado del eje X para evitar inestabilidad numérica ---
%===============================================================================
x_mean = mean(raman_shift);                    % Calcula el promedio del eje X
x_std = std(raman_shift);                      % Calcula la desviación estándar del eje X
x_scaled = (raman_shift - x_mean) / x_std;     % Normaliza los valores del eje X
%-------------------------------------------------------------------------------
%x_scaled = (raman_shift - x_mean) / x_std: Centra y escala los datos: 
%                resta el promedio → todos los datos giran en torno a 0.
%            Divide entre la desviación → todos los datos tienen unidad de dispersión.
%----------------------------------------------------------------------------------
%===============================================================================
% Ajuste polinomial de orden n=0,1,2,3... en datos escalados
%===============================================================================
orden_pol = 2;                                      % Define el orden del polinomio (manipula 0,1,2,3,4,5,6...)
p_scaled = polyfit(x_scaled, intensity, orden_pol); % Ajusta un polinomio a los datos
%------------------------------------------------------------------------------
   % x_scaled: eje X (datos escalados)
   % intensity: eje Y (intensidad del espectro)
   % orden_pol: orden del polinomio 
%-------------------------------------------------------------------------------
baselinea = polyval(p_scaled, x_scaled);             % Evalúa el polinomio (línea base)
   %p_scaled: polinomio que obtuviste con polyfit
   %x_scaled: valores de entrada (escalados)
%-------------------------------------------------------------------------------
%========================================================================
% Corrección del espectro
corrected_intensity = intensity - baselinea;
%-----------------------------------------------------------------
% --- Ajustar para que el mínimo sea cero ( Desplaza todo el espectro para que el valor más bajo sea 0) ---
corrected_intensity = corrected_intensity - min(corrected_intensity);
%-=======================================================================
% --- Graficar resultados -------------------------------------------
%-=======================================================================
figure;
plot(raman_shift, intensity, 'k', 'DisplayName','ZCe3.txt'); hold on; 
% Grafica el espectro original en negro
plot(raman_shift, baselinea, 'b--', 'LineWidth', 1.5, 'DisplayName','Línea Base');
% Grafica la línea base en azul punteado
plot(raman_shift, corrected_intensity, 'r', 'LineWidth', 1.2, 'DisplayName','Corregido (Linea Base)');
% Grafica el espectro corregido en rojo
%---------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
xlabel('Desplazamiento Raman (cm^{-1})');  % Etiqueta eje X
ylabel('Intensidad (u.a.)');               % Etiqueta eje Y
%title('Corrección de Línea Base - espectro Raman ZCe3'); % Título de la gráfica
legend('show'); grid on;                   % Muestra la leyenda y activa la cuadrícula
%================================================================================================
%================================================================================================
%%Guardar datos corregidos 
%================================================================================================
%================================================================================================
% --- Crear matriz de datos corregidos ---
datos_corregidos = [raman_shift, corrected_intensity];

% --- Ruta personalizada donde se guardará el archivo ---
carpeta_salida = 'C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Linea base';  % Cambia esta ruta si quieres otra

% Crear la carpeta si no existe
if ~exist(carpeta_salida, 'dir')
   mkdir(carpeta_salida);
end

% Nombre completo del archivo con la ruta
archivo_salida = fullfile(carpeta_salida, 'espectro_corregido.txt');

% Guardar archivo .txt separado por tabulaciones
writematrix(datos_corregidos, archivo_salida, 'Delimiter', '\t');

% Mensaje de confirmación
disp(['Espectro corregido guardado en: ', archivo_salida]);
%%%%%%%
exportgraphics(gcf, 'C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Resultados\Espectro_ZnO.png', 'Resolution', 300);
