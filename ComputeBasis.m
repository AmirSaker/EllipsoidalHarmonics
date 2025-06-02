function [D, constants] = ComputeBasis(v, a1, a2, a3, N)
N = 1;
num_points = size(v, 1);
num_harmonics = (N + 1)^2;
D = zeros(num_points, num_harmonics);
constants = zeros(num_harmonics, 1);

% Focal parameters
h1 = sqrt(a2^2 - a3^2);
h2 = sqrt(a1^2 - a3^2);
h3 = sqrt(a1^2 - a2^2);
a = [a1, a2, a3];

% Convert to ellipsoidal coordinates
[~, mu, nu] = Cart2Ellip(a, v);

col = 1;
    for n = 0:N
    % Get all Lamé functions for current degree
    [K_coeffs, ~] = computeLameK(h1, h2, h3, n);
    [L_coeffs, ~] = computeLameL(h1, h2, h3, n);
    [M_coeffs, ~] = computeLameM(h1, h2, h3, n);
    [N_coeffs, ~] = computeLameN(h1, h2, h3, n);
    
    % Combine all coefficients and classes
    all_coeffs = [K_coeffs; L_coeffs; M_coeffs; N_coeffs];
    all_classes = [repmat({'K'}, numel(K_coeffs), 1);
                 repmat({'L'}, numel(L_coeffs), 1);
                 repmat({'M'}, numel(M_coeffs), 1);
                 repmat({'N'}, numel(N_coeffs), 1)];

        for m = 1:2*n+1
        % Get current coefficients and class
        coeff = all_coeffs{m};
        class = all_classes{m};
        
        % Evaluate at rho = a1 (constant)
        E_rho = evalLame(coeff, class, a1, h2, h3);
        constants(col) = E_rho;
        
        % Vectorized evaluation for all points
        E_mu = arrayfun(@(x) evalLame(coeff, class, x, h2, h3), mu);
        E_nu = arrayfun(@(x) evalLame(coeff, class, x, h2, h3), nu);
        D(:, col) = E_rho .* E_mu .* E_nu;
        col = col + 1;
        end
    end
end