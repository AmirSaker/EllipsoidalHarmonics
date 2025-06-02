function value = evalLame(coeffs, class, x, h2, h3)
    % Evaluates Lamé functions at point x for a given class
    %
    % Inputs:
    %   coeffs : Polynomial coefficients [a_0, a_1, ..., a_n] (from computeLameX)
    %   class  : Lamé class ('K', 'L', 'M', or 'N')
    %   x      : Evaluation point(s) (can be vector)
    %   h2, h3 : Focal parameters (h? = ?(a?² - a?²), h? = ?(a?² - a?²))
    %
    % Output:
    %   value  : Evaluated Lamé function value(s)

    % Polynomial part evaluation (using Horner's method)
    poly_val = polyval(flip(coeffs), x); % polyval expects coefficients in descending order
        
    
    % Class-specific modifications
    switch upper(class)
        case 'K'
            % Pure polynomial: K(x) = P(x)
            value = poly_val;
            
        case 'L'
            % L(x) = sqrt|x² - h3²| = P(x)
            value = poly_val .* sqrt(abs(x.^2 - h3^2));
            
        case 'M'
            % M(x) = sqrt|x² - h2²| = P(x)
            value = poly_val .* sqrt(abs(x.^2 - h2^2));
            
        case 'N'
            % N(x) = sqrt|x² - h²|*sqrt|x² - h2²| =· P(x)
            value = poly_val .* sqrt(abs(x.^2 - h3^2)) .* sqrt(abs(x.^2 - h2^2));
            
        otherwise
            error('Unknown Lamé class. Use ''K'', ''L'', ''M'', or ''N''.');
    end
end