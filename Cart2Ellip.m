function [rho, mu, nu] = Cart2Ellip(a, x)
% Calcul des coefficients cubiques
    c2 = a(1).^2 + a(2).^2 + a(3).^2 - x(:,1).^2 - x(:,2).^2 - x(:,3).^2;
    c1 = a(1).^2*a(2).^2 + a(1).^2*a(3).^2 + a(2).^2*a(3).^2 ...
       - (a(2).^2 + a(3).^2)*x(:,1).^2 - (a(1).^2 + a(3).^2)*x(:,2).^2 - (a(1).^2 + a(2).^2)*x(:,3).^2;
    c0 = a(1).^2*a(2).^2*a(3).^2 - a(2).^2*a(3).^2*x(:,1).^2 - a(1).^2*a(3).^2*x(:,2).^2 - a(1).^2*a(2).^2*x(:,3).^2;
    
    % Résolution de l'équation cubique
    p = (c2.^2 - 3*c1) / 9;
    q = (9*c1.*c2 - 27*c0 - 2*c2.^3) / 54;
    s1=[];s2=[];s3=[];
    
    for i=1:length(p)
            if abs(q(i)) < abs(p(i)^(3/2))
        omega = acos(q(i)/p(i)^(3/2));
        else
           omega = 0;
        end
    
    % Calcul des racines
    s1(i) = 2*sqrt(p(i))*cos(omega/3) - c2(i)/3;
    s2(i) = 2*sqrt(p(i))*cos((omega - 2*pi)/3) - c2(i)/3;
    s3(i) = 2*sqrt(p(i))*cos((omega - 4*pi)/3) - c2(i)/3;
    
    % Trier chaque triplet (s1,s2,s3) par ligne
    all_roots = cat(2, s1(:), s2(:), s3(:));
    sorted_roots = sort(all_roots, 2, 'descend');
    
    % Conversion en coordonnées ellipsoïdales
%     rho(i) = sqrt(round(a(1).^2 + s1(i), 10));
%     mu (i)= sqrt(round(a(1).^2 + s2(i), 10));
%     nu (i)= sign(x(i,1))*sqrt(round(a(1).^2 + s3(i), 10));
%     
%     signs = sign(x);
    end
     rho = sqrt(round(a(1).^2 + sorted_roots(:,1), 10));
     mu = sqrt(round(a(1).^2 + sorted_roots(:,2), 10));
     nu = sign(x(i,1))*sqrt(round(a(1).^2 + sorted_roots(:,3), 10));
    
end