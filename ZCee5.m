% Definir la carpeta donde están los archivos de espectros 
carpeta = 'C:\Users\OVER REGINO\Desktop\Espectroscopía RAMAN - copia\Resultados';  % Carpeta donde se ubica el documento

% Obtener la lista de archivos .txt en la carpeta
archivos = dir(fullfile(carpeta, '*ZCe5.txt'));

% Verificar si hay archivos en la carpeta
if isempty(archivos)
    error('No se encontraron archivos .txt en la carpeta especificada.');
end

% Crear figura para visualizar todos los espectros
figure; hold on;
%title('Espectros Raman para el ZCe5'); 
xlabel('Desplazamiento Raman (cm^{-1})'); 
ylabel('Intensidad (u.a.)');

% Iterar sobre todos los archivos y graficarlos
for i = 1:length(archivos)
    archivo = fullfile(carpeta, archivos(i).name);
    
    % Leer datos desde el archivo de texto
    datos = readmatrix(archivo);

    % Verificar si el archivo tiene al menos dos columnas
    if size(datos, 2) < 2
        warning('El archivo %s no tiene el formato esperado. Se omitirá.', archivos(i).name);
        continue;
    end

    % Extraer columnas de desplazamiento Raman e intensidad
    desplazamientoRaman = datos(:,1);
    intensidad = datos(:,2);

    % Graficar cada espectro en negro con línea más gruesa
    plot(desplazamientoRaman, intensidad, 'k', 'LineWidth', 0.9, 'DisplayName', archivos(i).name);
end

% Mostrar leyenda con nombres de archivos
legend('show');
grid on;
%
exportgraphics(gcf, 'Espectro_ZnO.png', 'Resolution', 300);  % 300 dpi es calidad para publicación
