function [X,Q] = workspace()
    % butées mécaniques
    qmin =[-pi*120/180 -pi*120/180 -pi*180/180 -pi*60/180 -pi*90/180 -pi*70/180 -pi*20/180];
    qmax=[ pi*30/180 pi*50/180 pi*50/180 pi*90/180 pi*90/180 pi*90/180 pi*30/180];
    npoints = 15;  % on peut mettre plus de points pour un résultat plus précis.
    
    q_vectors = cell(4,1);
    for i = 1:4
        q_vectors{i} = linspace(qmin(i), qmax(i), npoints);
    end

    
    [q1,q2,q3,q4] = ndgrid(q_vectors{:});
    
    
    N = numel(q1);

    
    X = zeros(3, N);
    Q = zeros(7,N);
    
    for k = 1:N
        q = [q1(k); q2(k); q3(k); q4(k) ; 0 ; 0 ; 0];
        a = MGD_angle(q);
        X(:,k) = a(1:3);  % doit renvoyer [x;y;z]
        Q(:,k) = q;
    end
end

[X,Q] = workspace();  

figure;
scatter3(X(1,:), X(2,:), X(3,:), 10, 'filled');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
title('Espace de travail du robot');
axis equal;



% tracé coupe
x = X(1,:)';
y = X(2,:)';
z = X(3,:)';



z0 = 0.1;           % hauteur de coupe
dz = 0.02;          % tolérance autour de z0


idx = abs(X(3,:) - z0) < dz;


figure;
scatter(X(1,idx), X(2,idx), 20, 'filled');
xlabel('X'); ylabel('Y');
title(['Coupe horizontale à z = ', num2str(z0)]);
axis equal;
grid on;






% affichage manipulabilité
N = size(X, 2);
manip = zeros(1, N);

for k = 1:N
    qk = Q(:,k);
    J = jacobienne(qk);
     
    Jv = J(1:3, :); 
    
    manip(k) = sqrt(max(det(J*J'), 0));
end

figure('Color', 'w', 'Name', 'Analyse de Manipulabilité');
hold on;

scatter3(X(1,:), X(2,:), X(3,:), 15, manip, 'filled', 'MarkerFaceAlpha', 0.5);

shp = alphaShape(X(1,:)', X(2,:)', X(3,:)', 1.0);
h_shp = plot(shp, 'FaceColor', 'none', 'EdgeColor', [0.5 0.5 0.5], 'EdgeAlpha', 0.2);

colormap(jet); 
h_cb = colorbar;
ylabel(h_cb, 'Indice de Manipulabilité');
title('Gradient de Manipulabilité (Agilité du robot)');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

grid on; 
axis equal;
view(3);
camlight;



%tracé alphaShape (enveloppe lisse)
shp = alphaShape(x, y, z, 0.9);


figure;
plot(shp, 'FaceColor', 'cyan', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Enveloppe lisse de l''espace de travail (alpha shape)');
axis equal;
grid on;
view(3);
hLight = light('Position', [1 1 0], 'Style', 'infinite');

camlight headlight;      


rotate3d on;




