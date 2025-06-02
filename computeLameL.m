function [L_coeffs, p_values] = computeLameL(h1, h2, h3, n)
   % Calculate Lamé functions of class L
    % according to equations (22)-(24) from the PDF
    
    % Input validation
    if n < 0 || mod(n,1) ~= 0
        error('n must be a non-negative integer');
    end
    if h2 <= 0 || h3 <= 0
        error('h2 and h3 must be positive');
    end
    
    alpha = h3^2 + h2^2;
    beta = h3^2 * h2^2;
    
    % Determine matrix size based on n parity
    if mod(n, 2) == 0
        r = n/2;
    else
        r = (n+1)/2;
    end
    
    % Handle special case when r=0 (n=0)
    if r == 0
        L_coeffs = {1};
        p_values = [];
        return;
    end
    
    % Construct tridiagonal matrix L
    syms p;
    L_matrix = sym(zeros(r, r));
    
    for j = 1:r
        % Main diagonal
        L_matrix(j,j) = -alpha*(p - (n - 2*j + 1)^2) + (2*n - 4*j + 4)*h2^2;
        
        % Upper diagonal (j < r)
        if j < r
            L_matrix(j,j+1) = 2*j*(2*n - 2*j + 1);
        end
        
        % Lower diagonal (j > 1)
        if j > 1
            L_matrix(j,j-1) = -beta*(n - 2*j + 2)*(n - 2*j + 1);  % Fixed index
        end
    end
    
    % Solve det(L_matrix) = 0 for p
    try
        p_solutions = solve(det(L_matrix) == 0, p, 'Real', true);
        p_values = double(p_solutions);
        
        % Sort p_values for consistent output
        p_values = sort(p_values);
    catch
        warning('Failed to solve determinant equation');
        p_values = [];
        L_coeffs = {};
        return;
    end
    
    % Calculate coefficients for each p
    L_coeffs = cell(length(p_values), 1);
    
    for i = 1:length(p_values)
        current_p = p_values(i);
        b = zeros(r, 1);
        b(1) = 1; % Initial condition b0 = 1
        
        % Recursive coefficient calculation
        for k = 0:r-2
            numerator = (alpha*(current_p - (n - 2*k - 1)^2) - (2*n - 4*k - 1)*h2^2)*b(k+1);
            if k > 0
                numerator = numerator + beta*(n - 2*k + 2)*(n - 2*k + 1)*b(k);  % Fixed index
            end
            denominator = 2*(k + 1)*(2*n - 2*k - 1);
            
            if denominator == 0
                error('Division by zero in coefficient calculation');
            end
            b(k+2) = numerator / denominator;
        end
        
        L_coeffs{i} = b;
    end
end