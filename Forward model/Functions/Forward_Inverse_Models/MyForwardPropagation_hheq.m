function [Holo_est, uk]=...
    MyForwardPropagation_hheq(f,E,Nx,Ny,Nz,phase3D_H,phase3D_G,...
    born_order,F_nopad,F_pad,Ft,is_rm_alias,...
    Pupil,Avg_amp,isE2,isMemOpt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute internal 3D field using SEAGLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxiters = 10;
alpha = 0.005;

uk_minus1 = E;
uk_minus2 = E;
tk_minus1 = 1;

for k = 1:maxiters
   tk = 0.5 * (1+sqrt(1+tk_minus1^2));
   sk = uk_minus1 + ((tk_minus1 - 1)/tk) * (uk_minus1 - uk_minus2);
   [grad, b] = Grad_seagle(sk,E,f,phase3D_G,Nx,Ny,Nz,F_pad,Ft);
   uk = sk - alpha * grad;
   
   tk_minus1 = tk;
   uk_minus1 = uk;
   uk_minus2 = uk_minus1;
   
   cost = 0.5*norm(b(:), 2)^2;
   fprintf('[cost: %1.5f]\n',cost);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagate 3D field to image plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cEsp=sum(F(f.*E_k).*phase3D_H,3).*Pupil;
tmp = Ft(cEsp);
Holo_Efield_2D = tmp(1:Nx,1:Ny);

if isE2
    Holo_est = abs(Avg_amp+Holo_Efield_2D).^2;
else
    Holo_est = real(Holo_Efield_2D);
end
