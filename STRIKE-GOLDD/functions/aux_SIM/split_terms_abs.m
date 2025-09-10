function terms_abs = split_terms_abs(eq)
    % Convert the symbolic equation to a string and remove spaces
    equation_str = strrep(char(eq), ' ', '');
    
    % Initialize variables
    terms = {};
    start_idx = 1;
    parentheses = 0;
    length_eq = length(equation_str);
    
    % Process each character
    for i = 1:length_eq
        % Count parentheses to ignore signs inside them
        if equation_str(i) == '('
            parentheses = parentheses + 1;
        elseif equation_str(i) == ')'
            parentheses = parentheses - 1;
        end
        
        % Detect '+' or '-' only if we are not inside parentheses
        if parentheses == 0 && i > 1 && (equation_str(i) == '+' || equation_str(i) == '-')
            term = equation_str(start_idx:i-1);
            terms{end+1} = term;
            start_idx = i;
        end
    end
    
    % Add the last summand
    if start_idx <= length_eq
        term = equation_str(start_idx:end);
        terms{end+1} = term;
    end
    
    % Remove leading '+' or '-' signs to obtain absolute values
    terms_abs = cell(size(terms));
    for i = 1:length(terms)
        term = terms{i};
        if ~isempty(term) && (term(1) == '+' || term(1) == '-')
            term = term(2:end);
        end
        terms_abs{i} = term;
    end
end
