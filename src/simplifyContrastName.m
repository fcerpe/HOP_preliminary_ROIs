function contrastNameSimple = simplifyContrastName(contrastName)
    % Simplifies the contrast name by removing non-alphabetic characters and parentheses.
    %
    % Parameters:
    %   contrastName: The original contrast name as a string.
    %
    % Returns:
    %   contrastNameSimple: The simplified contrast name, with only lowercase alphabetic characters.
    %
    % Example usage:
    %   contrastNameSimple = simplifyContrastName('Contrast 1 (Session 1)');
    
    contrastNameSimple = regexprep(lower(contrastName), '[^a-z]+|\([^)]*\)', '');
end