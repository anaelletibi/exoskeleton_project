
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

    % --- Normaliser u et v ---
    u = u / norm(u);
    v = v / norm(v);

    % Vérifier orthogonalité (optionnel)
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

%________________________________________________________



function [h,h2]=bouger_robot_animer(q,q2,pas_affichage,temps)
    % Calcule et trace la trajectoire depuis les configurations articulaires
    
    % Supprimer q_init si présent
    

    
    % Échelles identiques sur X, Y et Z

    n = size(q,1);
  
    x=zeros(n-2,1);
    y=zeros(n-2,1);
    z=zeros(n-2,1);
    x2=zeros(n-2,1);
    y2=zeros(n-2,1);
    z2=zeros(n-2,1);
    
    for i = 1:n-2
        q_inter =  q(i+1,:);
        q_inter2 =  q2(i+1,:);
        
        matrices_inter = MGD_ortho(q_inter);
        x(i)=matrices_inter(1,4,7);
        y(i)=matrices_inter(2,4,7);
        z(i)=matrices_inter(3,4,7);
        matrices_inter2 = MGD_ortho(q_inter2);
        x2(i)=matrices_inter2(1,4,7);
        y2(i)=matrices_inter2(2,4,7);
        z2(i)=matrices_inter2(3,4,7);
    end
    %plot3(x-0.2, y, z, 'b-', 'LineWidth', 2, "HandleVisibility","off",'DisplayName', 'Trajectoire');

    %plot3(-x2-0.2, -0.3-y2, z2, 'b-', 'LineWidth', 2, "HandleVisibility","off",'DisplayName', 'Trajectoire');
    q_inter =  q(1,:);
    
    matrices_inter = MGD_ortho(q_inter);
    x7=squeeze(matrices_inter(1:3,1,7))/24;
    z7=squeeze(matrices_inter(1:3,3,7))/24;
        pts = [ [0;0;0] , squeeze(matrices_inter(1:3, 4, 4:5)),...
            squeeze(matrices_inter(1:3, 4, 5))+z7,...
            squeeze(matrices_inter(1:3, 4, 5))+z7/2,...
            squeeze(matrices_inter(1:3, 4, 5))+z7/2+x7,...
            squeeze(matrices_inter(1:3, 4, 5))+x7,...
            squeeze(matrices_inter(1:3, 4, 5)),...
            ];
    q2_inter =  q2(1,:);
    
    matrices_inter2 = MGD_ortho(q2_inter);
    x72=squeeze(matrices_inter2(1:3,1,7))/24;
    z72=squeeze(matrices_inter2(1:3,3,7))/24;
        pts2 = [ [0;0;0] , squeeze(matrices_inter2(1:3, 4, 4:5)),...
            squeeze(matrices_inter2(1:3, 4, 5))+z72,...
            squeeze(matrices_inter2(1:3, 4, 5))+z72/2,...
            squeeze(matrices_inter2(1:3, 4, 5))+z72/2+x72,...
            squeeze(matrices_inter2(1:3, 4, 5))+x72,...
            squeeze(matrices_inter2(1:3, 4, 5)),...
            ];
    h=plot3(pts(1,:)-0.2, pts(2,:)+0.3, pts(3,:), '-o', 'Color', [0 0 0], ...
          'LineWidth', 2, 'MarkerSize', 4, 'HandleVisibility', 'off');
    h2=plot3(-pts2(1,:)-0.2, -pts2(2,:)-0.3, pts2(3,:), '-o', 'Color', [0 0 0], ...
          'LineWidth', 2, 'MarkerSize', 4, 'HandleVisibility', 'off');
    
      for i = 2:pas_affichage:n-1
            q_inter =  q(i,:);
            q_inter(7)=q_inter(7) - pi/2;
            
            matrices_inter = MGD_ortho(q_inter);
            x7=squeeze(matrices_inter(1:3,1,7))/24;
        z7=squeeze(matrices_inter(1:3,3,7))/24;
        pts = [ [0;0;0] , squeeze(matrices_inter(1:3, 4, 4:5)),...
            squeeze(matrices_inter(1:3, 4, 5))+z7,...
            squeeze(matrices_inter(1:3, 4, 5))+z7/2,...
            squeeze(matrices_inter(1:3, 4, 5))+z7/2+x7,...
            squeeze(matrices_inter(1:3, 4, 5))+x7,...
            squeeze(matrices_inter(1:3, 4, 5)),...
            ];
    q2_inter =  q2(i,:);
            q2_inter(7)=q2_inter(7) - pi/2;

    
    matrices_inter2 = MGD_ortho(q2_inter);
    x72=squeeze(matrices_inter2(1:3,1,7))/24;
    z72=squeeze(matrices_inter2(1:3,3,7))/24;
        pts2 = [ [0;0;0] , squeeze(matrices_inter2(1:3, 4, 4:5)),...
            squeeze(matrices_inter2(1:3, 4, 5))+z72,...
            squeeze(matrices_inter2(1:3, 4, 5))+z72/2,...
            squeeze(matrices_inter2(1:3, 4, 5))+z72/2+x72,...
            squeeze(matrices_inter2(1:3, 4, 5))+x72,...
            squeeze(matrices_inter2(1:3, 4, 5)),...
            ];
            set(h, 'XData', pts(1,:)-0.2);
            set(h, 'YData', pts(2,:)+0.3);
            set(h, 'ZData', pts(3,:));
            set(h2, 'XData', -pts2(1,:)-0.2);
            set(h2, 'YData', -pts2(2,:)-0.3);
            set(h2, 'ZData', pts2(3,:));
            pause(temps)
      end

    end


T = 10;
temps = linspace(0, T, 100);
% Définir les points de départ et d'arrivée


% Définir les points de départ et d'arrivée


p0 = [0,0,-0.5,pi/4,pi/4,0;...
    0,-0.3,-0.3,0,0, -pi/4;...
    0,-0.3,-0.3, 0,-pi/4,0 ;...
    0,-0.1,-0.3,0, pi/4,0;...
    0,-0.1,-0.3,0,pi/4,0;...
    0,-0.1,-0.3,0, pi/4,0];
pf = [-0.2,0,-0.4, 0,0,-pi/4;...
    -0.2,0,-0.4, 0,-pi/4,0;...
    0,-0.1,-0.3,0, pi/4,0;...
    0,-0.25,-0.3,0, pi/4,0;...
    0,-0.25,-0.3,0,pi/4,0;...
    0,0,-0.5,pi/4,pi/4,0;];
type_courbe=[0;1;0;0;-1;0];
b=[0;0.2;0;0;0.05;0];

p02 = [0,0,-0.5,pi/4,pi/4,0;...
    0,-0.3,-0.3,0,0, -pi/4;...
    0,-0.3,-0.3,0, pi/4,0;...
    0,-0.1,-0.3,0,pi/4,0;...
    0,-0.1,-0.3, 0,pi/4,0;...
    0,-0.1,-0.3,0,pi/4,0];
pf2 = [-0.2,0,-0.4, 0,0,-pi/4;...
    -0.2,0,-0.4, 0,pi/4,0;...
    0,-0.1,-0.3,0,pi/4,0;...
    0,-0.25,-0.3,0, pi/4,0;...
    0,-0.25,-0.3,0,pi/4,0;...
    0,0,-0.5,pi/4,pi/4,0];
type_courbe2=[0;1;0;0;-1;0];
b2=[0;0.2;0;0;0.05;0];

figure('Position', [100, 100, 800, 800]);
hold on

N=size(p0,1);
%[p_desired, v_desired] = quasi_rectiligne_cloche(p0, pf, temps, T);

        q_init = MGI_NewtonRaphson_contraintes(p0(1,:),[ 0 0 0 0 0 0 0]);
    for i=1:N
        if abs(type_courbe(i))==1
            
            u = [pf(i,1)-p0(i,1),pf(i,2)-p0(i,2),pf(i,3)-p0(i,3)];
            v = [-u(2),u(1),0];
            sens = type_courbe(i);
            [p_desired, v_desired] = quasi_ellipse_cloche(v,b(i),p0(i,:),pf(i,:), sens, temps, T);
        else
            [p_desired, v_desired] = quasi_rectiligne_cloche(p0(i,:), pf(i,:), temps, T);
        
        end

        q = q_depuis_trajectoire(p_desired, q_init);
        plot(linspace(-99+100*i,+100*i,101),q)
        q_init=q(end,:);
    end
hold off

figure('Position', [100, 100, 800, 800]);
    
    view(3); xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

hold on
xlim([-0.6 0.5])
ylim([-0.6 0.5])
zlim([-0.6 0.5])

xlabel('X(m)');
ylabel('Y(m)');
zlabel('Z(m)');
title('trajectoire de "check"généré pour deux exosquelette');
legend('show');
grid on;
q_init = MGI_NewtonRaphson_contraintes(p0(1,:),[ 0 0 0 0 0 0 0]);
q_init2 = MGI_NewtonRaphson_contraintes(p02(1,:),[ 0 0 0 0 0 0 0]);
pause(3)

for j=1:10
    for i=1:N
        q_init = MGI_NewtonRaphson_contraintes(p0(i,:),[ 0 0 0 0 0 0 0]);
        if abs(type_courbe(i))==1

            u = [pf(i,1)-p0(i,1),pf(i,2)-p0(i,2),pf(i,3)-p0(i,3)];
            v = [-u(2),u(1),0];
            sens = type_courbe(i);
            [p_desired, v_desired] = quasi_ellipse_cloche(v,b(i),p0(i,:),pf(i,:), sens, temps, T);
        else
            [p_desired, v_desired] = quasi_rectiligne_cloche(p0(i,:), pf(i,:), temps, T);
        
        end
        q = q_depuis_trajectoire(p_desired, q_init);
        
        
        if abs(type_courbe2(i))==1
            
            u = [pf2(i,1)-p02(i,1),pf2(i,2)-p02(i,2),pf2(i,3)-p02(i,3)];
            v = [-u(2),u(1),0];
            sens = type_courbe2(i);
            [p_desired, v_desired] = quasi_ellipse_cloche(v,b2(i),p02(i,:),pf2(i,:), sens, temps, T);
        else
            [p_desired, v_desired] = quasi_rectiligne_cloche(p02(i,:), pf2(i,:), temps, T);
        
        end
        q2 = q_depuis_trajectoire(p_desired, q_init2);
        [h,h2]=bouger_robot_animer(q(:,1:7),q2(:,1:7),3,0.08);
        delete(h)
        delete(h2)
        q_init=q(end,:);
        q_init2=q2(end,:);
    end
end
hold off