function [class, idx] = LameIndex(n, m)
    % Détermine la classe (K,L,M,N) et l'index pour m donné
    % Input:
    %   n : degré de la fonction de Lamé (0, 1, 2,...)
    %   m : ordre (1 à 2n+1)
    % Output:
    %   class : 1 (K), 2 (L), 3 (M), 4 (N)
    %   idx : index dans la classe
    
    % Nombre de fonctions par classe pour le degré n
    class_sizes = [n - 1,  % K
                   n-2,      % L
                   n-2,      % M
                   max(n - 3, 0)]; % N
    
    % Détermination de la classe
    if m <= sum(class_sizes(1))  % Classe K
        class = 1;
        idx = m;
    elseif m <= sum(class_sizes(1:2))  % Classe L
        class = 2;
        idx = m ;
    elseif m <= sum(class_sizes(1:3))  % Classe M
        class = 3;
        idx = m ;
    else  % Classe N
        class = 4;
        idx = m ;
    end
end