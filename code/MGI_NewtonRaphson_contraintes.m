function q = MGI_NewtonRaphson_contraintes(xd , q0)

    maxIter = 500;
    eps = 1e-4;
    
    
    q = q0(:);
    xd = xd(:);
    % butées mécaniques
    qmin =[-pi*120/180 -pi*120/180 -pi*180/180 -pi*60/180 -pi*90/180 -pi*70/180 -pi*20/180];
    qmax=[ pi*30/180 pi*50/180 pi*50/180 pi*90/180 pi*90/180 pi*90/180 pi*30/180];
    qmean = (qmax + qmin)/2;
    
    
    for k = 1:maxIter

        T = MGD_ortho(q);
        T07 = T(:,:,7);

        R = T07(1:3,1:3);
        p = T07(1:3,4);

        phi   = atan2(R(3,2), R(3,3));
        theta = -asin(R(3,1));
        psi   = atan2(R(2,1), R(1,1));

        x = [p; phi; theta; psi];

        e = xd - x;
        
        e(4:6) = wrapToPi(e(4:6));
        if norm(e) < eps  % on sort de la boucle si on a une précision suffisante
            fprintf('Convergence atteinte en %d iterations\n', k);
            return
        end

        J = my_exo_orth_jac(q);
        Jp = pinv(J);

        delta_q0 = - (q - qmean) ./ ( (qmax - qmin).^2 );
        m =eye(7)-Jp*J;
        dq = 0.5*(Jp*e + (m)*delta_q0) ;
        q = q + dq;
        q = q(:, end);
        
        

    end

    warning('MGI : convergence non atteinte');
    
end
