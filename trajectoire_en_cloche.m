 function tracer_position_3d(p_trajectory)
    % Trace la position p en 3 dimensions
    % Args:
    %   p_trajectory: Matrice contenant les positions [x, y, z] au fil du temps
    figure('Position', [100, 100, 800, 800]);
    
    if iscell(p_trajectory)
        p_trajectory = cell2mat(p_trajectory');
    end
    
    x = p_trajectory(:, 1);
    y = p_trajectory(:, 2);
    z = p_trajectory(:, 3);
    
    plot3(x, y, z, 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectoire générée');
    hold on;
    scatter3(x(1), y(1), z(1), 100, 'green', 'filled', 'DisplayName', 'Position initiale');
    scatter3(x(end), y(end), z(end), 100, 'red', 'filled', 'DisplayName', 'Position finale');
    
    % % trajectoire reciligne idéale
    line([x(1) x(end)], [y(1) y(end)], [z(1) z(end)], ...
        'Color', 'w', 'LineStyle', '--', 'LineWidth', 2, ...
        'DisplayName', 'Trajectoire en commande')

    figure;
    plot3(x, y, z,'w', 'LineWidth', 2);
    grid on;

    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    title('Trajectoire quasi-rectiligne');
    legend('show');
    grid on;
    
    
    all_coords = [x; y; z];
    coord_min = min(all_coords);
    coord_max = max(all_coords);
    
    
    marge = 0.2;  % 20% de marge
    range = coord_max - coord_min;
    coord_min = coord_min - range * marge;
    coord_max = coord_max + range * marge;
    
   
    xlim([coord_min, coord_max]);
    ylim([coord_min, coord_max]);
    zlim([coord_min, coord_max]);
    
    
    axis equal;
    
    hold off;
end

function dx = derivee(x, t)
    
    dx = zeros(length(x), 1);
    for i = 1:length(x)-1
        dx(i) = (x(i+1) - x(i)) / (t(i+1) - t(i));
    end
end

function S = integrale(x, t)
    
    S = zeros(length(x), 1);
    S(1) = 0;
    for i = 1:length(x)-1
        S(i+1) = S(i) + (x(i) + x(i+1)) * (t(i+1) - t(i)) / 2;
    end
end

function [p,v]= bouger(q, t)
    % Calcule et trace la trajectoire depuis les configurations articulaires
    if size(q,1) == length(t) + 1
        q = q(2:end, :);
    end
    [T, n] = size(q);
    
    % Calcul des positions avec my_exo_orth_jac
    pos_x = zeros(T, 1);
    pos_y = zeros(T, 1);
    pos_z = zeros(T, 1);
    
    p = zeros(length(t), 3);

    for k = 1:T
        % Calcul des matrices homogènes pour la configuration q(k,:)
        T_matrices = MGD_ortho(q(k,:)');
        
        % Extraction de la position de l'effecteur terminal (T07)
        T07 = T_matrices(:,:,7);
        pos_x(k) = T07(1,4); 
        pos_y(k) = T07(2,4);
        pos_z(k) = T07(3,4);
        p(k,1) = pos_x(k);
        p(k,2) = pos_y(k);
        p(k,3) = pos_z(k);
    end
    
    % Assemblage de la trajectoire
    p_trajectory = [pos_x, pos_y, pos_z];

        % Tracé de la trajectoire
    tracer_position_3d(p_trajectory);
    
    % Calcul des vitesses par dérivation 
    v_x = derivee(pos_x, t);
    v_y = derivee(pos_y, t);
    v_z = derivee(pos_z, t);
    v = [v_x, v_y, v_z];
    
    end



function [x,v] = quasi_rectiligne_cloche(p0, pf, temps, T)
    x = zeros(length(temps), 6);
    v = zeros(length(temps), 6);
    
    for i = 1:length(temps)
        t = temps(i);
        s = t/T;

        profil = 10*s^3 - 15*s^4 + 6*s^5;
        dprofil = (30*s^2 - 60*s^3 + 30*s^4)/T;

        x(i, 1) = p0(1) + (pf(1) - p0(1)) * profil;
        x(i, 2) = p0(2) + (pf(2) - p0(2)) * profil;
        x(i, 3) = p0(3) + (pf(3) - p0(3)) * profil;
        x(i, 4) = p0(4) + (pf(4) - p0(4)) * profil;
        x(i, 5) = p0(5) + (pf(5) - p0(5)) * profil;
        x(i, 6) = p0(6) + (pf(6) - p0(6)) * profil;

        v(i,1) = (pf(1) - p0(1)) * dprofil;
        v(i,2) = (pf(2) - p0(2)) * dprofil;
        v(i,3) = (pf(3) - p0(3)) * dprofil;
        v(i,4) = (pf(4) - p0(4)) * dprofil;
        v(i,5) = (pf(5) - p0(5)) * dprofil;
        v(i,6) = (pf(6) - p0(6)) * dprofil;
    end
end

function [x,v] = quasi_ellipse_cloche(v, b, p0, pf, sens, temps,T)

    % C : centre [Cx Cy Cz]
    % u : direction axe a (vecteur 3D)
    % v : direction axe b (vecteur 3D)
    % a : demi-grand axe
    % b : demi-petit axe
    % sens : +-1 pour sens trigonométrique/non trigonométrique

    C = [(p0(1)+pf(1))/2,(p0(2)+pf(2))/2,(p0(3)+pf(3))/2];
    a = sqrt((p0(1)-pf(1))^2+(p0(2)-pf(2))^2+(p0(3)-pf(3))^2)/2;

    u = [pf(1)-p0(1),pf(2)-p0(2),pf(3)-p0(3)];

    % Normaliser u et v
    u = u / norm(u);
    v = v / norm(v);

    % Vérifier orthogonalité
    if abs(dot(u,v)) > 1e-6
        error('u et v doivent être orthogonaux');
    end

    theta = zeros(length(temps), 1);
    dtheta = zeros(length(temps), 1);

    alpha = zeros(length(temps),3);
    dalpha = zeros(length(temps),3);

    for i = 1:length(temps)
        t = temps(i);
        s = t/T;
        profil = 10*s^3 - 15*s^4 + 6*s^5;
        dprofil = (30*s^2 - 60*s^3 + 30*s^4)/T;
        theta(i) = pi * profil;
        dtheta(i) = pi * dprofil;

        alpha(i, 1) = p0(4) + (pf(4) - p0(4)) * profil;
        alpha(i, 2) = p0(5) + (pf(5) - p0(5)) * profil;
        alpha(i, 3) = p0(6) + (pf(6) - p0(6)) * profil;

        dalpha(i,4) = (pf(4) - p0(4)) * dprofil;
        dalpha(i,5) = (pf(5) - p0(5)) * dprofil;
        dalpha(i,6) = (pf(6) - p0(6)) * dprofil;
    end

    % Paramétrisation de l'ellipse
    P = C + a*cos(sens*theta).*u + b*sin(sens*theta).*v;
    vP = -a*sens*dtheta.*sin(sens*theta).*u+b*sens*dtheta.*cos(sens*theta).*v;

    x=[P(:,1),P(:,2),P(:,3),alpha(:,1),alpha(:,2),alpha(:,3)];
    v=[vP(:,1),vP(:,2),vP(:,3),dalpha(:,1),dalpha(:,2),dalpha(:,3)];

    % Tracé
    plot3(x(:,1), x(:,2), x(:,3), 'LineWidth', 2)
    hold on
    scatter3(C(1), C(2), C(3), 80, 'filled')
    
    xlabel('X'); ylabel('Y'); zlabel('Z')
    grid on
    axis equal
end

function q_list = q_depuis_trajectoire(p, q_init)
    q = q_init;
    nSteps = size(p, 1);
    nq = numel(q_init);
    q_list = zeros(nSteps+1, nq);
    q_list(1, :) = q_init(:).';
    for i = 1:nSteps
        x = p(i, :);
        q = MGI_NewtonRaphson_contraintes(x, q);
        q_list(i+1, :) = q(:).';
    end
end



T = 10;
temps = linspace(0, T, 100);


p0 = [0.100, 0.150, 0.200, 0, 0, 0];
pf = [0.350, 0.200, 0.250, 0, 0, 0];
q_init = [0, 0, 0, 0, 0, 0, 0];




b=0.2;
u = [pf(1)-p0(1),pf(2)-p0(2),pf(3)-p0(3)];
v = [-u(2),u(1),0];
sens = 1;

[p_desired, v_desired] = quasi_rectiligne_cloche(p0, pf, temps, T);


q = q_depuis_trajectoire(p_desired, q_init);


[p,v] = bouger(q, temps);



%tracé de la vitesse et de la position en fonction du temps
figure
hold on
plot(temps, p_desired(:,1),'w','LineWidth',2)
plot(temps, p_desired(:,2), 'r','LineWidth',2)
plot(temps, p_desired(:,3), 'y','LineWidth',2)
plot(temps, p(:,1), 'x', 'MarkerSize', 5, 'MarkerEdgeColor', 'w')
plot(temps, p(:,2), 'x', 'MarkerSize', 5, 'MarkerEdgeColor', 'r')
plot(temps, p(:,3), 'x', 'MarkerSize', 5, 'MarkerEdgeColor', 'y')
hold off

grid on
xlabel('Temps (s)')
ylabel('Position (m)')
legend('x commande','y commande','z commande','x généré', 'y généré', 'z généré')
title('Comparaison des positions du bout de l effecteur en commande et générée')

figure
plot(temps, v_desired(:,1),'w', 'LineWidth',2)
hold on
plot(temps, v_desired(:,2),'r', 'LineWidth',2)
plot(temps, v_desired(:,3),'y', 'LineWidth',2)
plot(temps, v(:,1), 'x', 'MarkerSize', 5, 'MarkerEdgeColor', 'w')
plot(temps, v(:,2), 'x', 'MarkerSize', 5, 'MarkerEdgeColor', 'r')
plot(temps, v(:,3), 'x', 'MarkerSize', 5, 'MarkerEdgeColor', 'y')
grid on
xlabel('Temps (s)')
ylabel('Vitesse (m/s)')
legend('vx commande','vy commande','vz commande','vx généré', 'vy généré', 'vz généré')
title('Comparaison des vitesses du bout de l effecteur idéales et simulées pour une loi en cloche')


% erreur relative

% Différence position
dx = p(:,1) - p_desired(:,1);
dy = p(:,2) - p_desired(:,2);
dz = p(:,3) - p_desired(:,3);

% Norme erreur
erreur = 100*sqrt(dx.^2 + dy.^2 + dz.^2);

% Norme commande
norm_cmd = sqrt(p_desired(:,1).^2 + p_desired(:,2).^2 + p_desired(:,3).^2);

% Éviter division par zéro
norm_cmd(norm_cmd < 1e-8) = 1e-8;

% Erreur relative
erreur_relative = erreur ./ norm_cmd;

% Tracé
figure
plot(temps, erreur_relative, 'w', 'LineWidth', 2)
grid on
xlabel('Temps (s)')
ylabel('Erreur relative (en %) ')
title('Erreur relative entre la commande en position et la position générée')



% Différence vitesse
dx = v(:,1) - v_desired(:,1);
dy = v(:,2) - v_desired(:,2);
dz = v(:,3) - v_desired(:,3);

% Norme erreur
erreur = 100*sqrt(dx.^2 + dy.^2 + dz.^2);

% Norme commande
norm_cmd = sqrt(v_desired(:,1).^2 + v_desired(:,2).^2 + v_desired(:,3).^2);

% Éviter division par zéro
norm_cmd(norm_cmd < 1e-8) = 1e-8;

% Erreur relative
erreur_relative = erreur ./ norm_cmd;

% Tracé
figure
plot(temps, erreur_relative, 'w', 'LineWidth', 2)
grid on
xlabel('Temps (s)')
ylabel('Erreur relative (en %) ')
title('Erreur relative entre la commande en vitesse et la vitesse générée')
