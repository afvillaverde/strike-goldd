function sumandos = split_terms_abs(ecuacion)
    % separa los sumandos de una ecuación y los devuelve en un cell array
    
    % Eliminar espacios en blanco para facilitar el procesamiento
    ecuacion = strrep(char(ecuacion), ' ', '');
    
    % Inicializar variables
    sumandos = {};
    inicio = 1;
    longitud = length(ecuacion);
    
    % Procesar cada carácter de la ecuación
    for i = 1:longitud
        % Detectar cambio de signo (excepto en el primer carácter)
        if i > 1 && (ecuacion(i) == '+' || ecuacion(i) == '-')
            % Guardar el sumando anterior (con su signo)
            sumando = ecuacion(inicio:i-1);
            sumandos{end+1} = sumando;
            inicio = i+1;
        end
    end
    
    % Añadir el último sumando
    if inicio <= longitud
        sumando = ecuacion(inicio:end);
        sumandos{end+1} = sumando;
       
    end
    
% Eliminar signos para obtener valor absoluto
    sumandos_abs = cell(size(sumandos));
    for i = 1:length(sumandos)
        sumando = sumandos{i};
        % Eliminar signo si existe al inicio
        if sumando(1) == '+' || sumando(1) == '-'
            sumando = sumando(2:end);
        end
        sumandos_abs{i} = sumando;
    end
end