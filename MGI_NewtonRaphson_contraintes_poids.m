function q = MGI_NewtonRaphson_contraintes_poids(xd, q0)
    maxIter = 500;
    eps = 1e-6;
    q = q0(:);
    xd = xd(:);
    
    qmin = [-120, -120, -180, -60, -90, -70, -20]' * pi/180;
    qmax = [  30,   50,   50,  90,  90,  90,  30]' * pi/180;
    qmean = (qmax + qmin) / 2;

    
    weights = [0.2, 0.2, 0.2, 10, 10, 10, 10];
    W = diag(weights);
    Winv = inv(W);

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
        e(4:6) = (atan2(sin(e(4:6)), cos(e(4:6)))); 
        
        if norm(e) < eps
            fprintf('Convergence atteinte en %d iterations\n', k);
            return
        end
        
        J = my_exo_orth_jac(q);
        
        
        Jp_w = Winv * J' / (J * Winv * J');
        
        
        grad_confort = (q - qmean) ./ (qmax - qmin).^2;
        
        I = eye(7);
        dq = Jp_w * e - (I - Jp_w * J) * grad_confort * 0.1;
        
        q = q + dq;
    end
    warning('MGI : convergence non atteinte');
end