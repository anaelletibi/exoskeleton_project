% Calcule la jacobienne pour le MGD à axes orthogonaux

function J = my_exo_orth_jac(q)
    matrices = MGD_ortho(q);
    n = size(matrices,3);

    vecteur_end = matrices(1:3,4,n);   % OoOend

    Jv = zeros(3,n);   % initialisation de Jv
    Jw = zeros(3,n);     % initialisation de Jw

    for i = 1:n
        z_i = matrices(1:3,3,i);
        vecteur_i = matrices(1:3,4,i);

        Jv(:,i) = cross(z_i, vecteur_end - vecteur_i) ;     % on rajoute les dérivées
        Jw(:,i) = z_i;                       % on rajoute les zi
    end

    J = [Jv; Jw];                % on concatène
end
