function [Holo_est, uk]=...
    MyForwardPropagation_hheq_btrack(f,E,Nx,Ny,Nz,phase3D_H,phase3D_G,...
    F_pad,Ft,Pupil,Avg_amp,isE2,holo_b1,holo_bm,result_dir,maxiters,...
    save_every_n_iter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute internal 3D field using SEAGLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step = 0.1;
backTol = 1e-12;
eta = 1/2;
is_save_intermediate_data = 1;

uk_minus1 = E;
uk_minus2 = E;
tk_minus1 = 1;

cost_vec = zeros(1,maxiters);
for k = 1:maxiters
    tic
    tk = 0.5 * (1+sqrt(1+tk_minus1^2));
    sk = uk_minus1 + ((tk_minus1 - 1)/tk) * (uk_minus1 - uk_minus2);
    [grad, resid] = Grad_seagle(sk,E,f,phase3D_G,Nx,Ny,Nz,F_pad,Ft);
    uk = sk - step * grad;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All additional code for backtracking goes here (waleed)
    d = 0.5*norm(resid(:), 2)^2;
    uk_resid = Seagle_A_operator(uk,E,f,phase3D_G,Nx,Ny,Nz,F_pad,Ft) - E;
    data = 0.5*norm(uk_resid(:), 2)^2;
    quad = d...
        + (0.5/step)*(norm(uk(:)-sk(:)+step*grad(:))^2)...
        - (0.5*step)*(norm(grad(:))^2);
    
    % Backtracking check the line search condition
    while(backTol > 0 && (quad - data) < 0)
        % reduce step
        step = eta*step;
        
        % tentative update
        %xhatnext = proxTV(s-step*grad, step*tau, gradObj);
        uk = sk - step * grad;
        
        % data-term and quadratic upper bound
        %[data,z_est,SF,E_rb_all] = evaluateFidCost(xhatnext,Holo2D,forwardObj);
        uk_resid = Seagle_A_operator(uk,E,f,phase3D_G,Nx,Ny,Nz,F_pad,Ft) - E;
        data = 0.5*norm(uk_resid(:), 2)^2;
        quad = d...
            + (0.5/step)*(norm(uk(:)-sk(:)+step*grad(:))^2)...
            - (0.5*step)*(norm(grad(:))^2);
        
        % print to command line
        if(1)
            fprintf('[Iter: %d]', k);
            fprintf('[update attempt]');
            fprintf('[step: %.1e]', step);
            fprintf('[data: %.1e]', data);
            fprintf('[quad: %.1e]', quad);
            fprintf('[diff: %.1e]', quad-data);
            fprintf('\n');
        end
        
        % if step is smaller than tolerance
        if(step < backTol)
            warning(' step < %.1e', backTol);
            step = backTol;
            break;
        end
     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tk_minus1 = tk;
    uk_minus1 = uk;
    uk_minus2 = uk_minus1;
    
    cost_vec(k) = 0.5*norm(resid(:), 2)^2;
    fprintf('[iter: %d] [cost: %1.5f] [iter_time: %f]\n',k,cost_vec(k),...
        toc);
    
    if(~rem(k,save_every_n_iter) && is_save_intermediate_data)
        save_intermediate_data;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagate 3D field to image plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cEsp=sum(F_pad(f.*uk).*phase3D_H,3).*Pupil;
tmp = Ft(cEsp);
Holo_Efield_2D = tmp(1:Nx,1:Ny);

if isE2
    Holo_est = abs(Avg_amp+Holo_Efield_2D).^2;
else
    Holo_est = real(Holo_Efield_2D);
end
