function [D, constants] = computeSurfaceBasis(v, a1, a2, a3, N)

   % Initialize
    num_points = size(v, 1);
    num_harmonics = (N + 1)^2;
    D = zeros(num_points, num_harmonics);
    constants = zeros(num_harmonics, 1); % Stores E_n^m(a1)
    
    % Focal parameters
    h1 = sqrt(a2^2 - a3^2);
    h2 = sqrt(a1^2 - a3^2);
    h3 = sqrt(a1^2 - a2^2);
    a=[a1,a2,a3];
    x=[v(:,1), v(:,2), v(:,3)];
    % Convert to ellipsoidal coordinates (?,? only, ?=a1)
    [~, mu, nu] = Cart2Ellip(a,x);
    
    col = 1;
    for n = 0:N
        % Get Lamé functions for current degree
        [K_coeffs, ~] = computeLameK(h1, h2, h3, n);
        [L_coeffs, ~] = computeLameL(h1, h2, h3, n);
        [M_coeffs, ~] = computeLameM(h1, h2, h3, n);
        [N_coeffs, ~] = computeLameN(h1, h2, h3, n);
        
        %test
        coeff=cat(1,K_coeffs,L_coeffs,M_coeffs,N_coeffs);
        
        for m = 1:2*n+1
            [class, idx] = getLameIndex(n, m);
            % Evaluate ?-component (constant on surface)
            switch class
                case 1, E_rho = evalLame(coeff{idx,1}, 'K', a1, h2, h3);
                case 2, E_rho = evalLame(coeff{idx,1}, 'L', a1, h2, h3);
                case 3, E_rho = evalLame(coeff{idx,1}, 'M', a1, h2, h3);
                case 4, E_rho = evalLame(coeff{idx,1}, 'N', a1, h2, h3);
            end
            constants(col) = E_rho;
            
            % Evaluate and components
            for i = 1:num_points
              for m = 1:2*n+1
                [class, idx] = getLameIndex(n, m);
                switch class
                    case 1
                        E = E_rho * evalLame(coeff{idx,1}, 'K', mu(i), h2, h3) ...
                              * evalLame(coeff{idx,1}, 'K', nu(i), h2, h3);
                    case 2
                        E = E_rho * evalLame(coeff{idx,1}, 'L', mu(i), h2, h3) ...
                              * evalLame(coeff{idx,1}, 'L', nu(i), h2, h3);
                    case 3
                        E = E_rho * evalLame(coeff{idx,1}, 'M', mu(i), h2, h3) ...
                              * evalLame(coeff{idx,1}, 'M', nu(i), h2, h3);
                    case 4
                        E = E_rho * evalLame(coeff{idx,1}, 'N', mu(i), h2, h3) ...
                              * evalLame(coeff{idx,1}, 'N', nu(i), h2, h3);
                 end
                          
                end
                D(i, col) = E;
            end
            col = col + 1;
        end
    end

end