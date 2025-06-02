function [x, y, z] = ellipsoidal2Cart(a1, a2, a3, rho, mu, nu)
    % Convertit les coordonnées ellipsoïdales en cartésiennes

    % Excentricités linéaires
    h1 = sqrt(a2^2 - a3^2);
    h2 = sqrt(a1^2 - a3^2);
    h3 = sqrt(a1^2 - a2^2);

    % Coordonnées cartésiennes
    x = (rho * mu * nu) / (h2 * h3);
    y = sqrt(rho^2 - h3^2) * sqrt(mu^2 - h3^2) * sqrt(h3^2 - nu^2) / (h1 * h3);
    z = sqrt(rho^2 - h2^2) * sqrt(h2^2 - mu^2) * sqrt(h2^2 - nu^2) / (h1 * h2);
end
