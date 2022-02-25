function [Holo_est, Holo_Efield_2D, E_rb_all,Perf_all]=...
    MyForwardPropagation_rb3D(f,E,Nx,Ny,Nz,phase3D_H,phase3D_G,...
    born_order,F_nopad3D,F_pad3D,F2D,Fpad2D,Ft3D,Ft2D,is_rm_alias,...
    Pupil,Avg_amp,isE2,isMemOpt)

%%%
Nzg = 2*(Nz+1)-1; %numslices in Phase3D_G
Nzg_mid = (Nzg+1)/2; %idx of zero-propagation slice
%%%
%
if is_rm_alias
    cEsp=zeros(2*Nx,2*Ny,Nz);
    F = F_pad3D;
    F2 = Fpad2D;
else
    cEsp=zeros(Nx,Ny,Nz);
    F = F_nopad3D;
    F2 = Fpad2D;
end

% This matrix saves all born fields uk, uk-1, ... , u0
E_rb_all = zeros(Nx,Ny,Nz,born_order);
E_rb_all(:,:,:,1) = E;


%% Implementation of Eq 7b from Ulgubek Recursive Born
%  Recursive computation of kth order born field within the object
E_k = E;
E_s  = E;

k = 1;
while(k < born_order)
    
    % 1 Multiply k-1th order born scatteredfield with the object
    %Ef=E_s.*f;
    
    % Compute 2D FFTs for all slices of the previous product
    % Compute FFT of object-incident_field product
    % (for filtering in fourier domain)
    
    F_Ef = F(E_s.*f);
    cEsp = F_Ef.*phase3D_G;   
    tmp = Ft3D(cEsp);
    E_s = tmp(1:Nx,1:Ny,1:Nz);
    clear tmp;
    
    E_k = E_k + E_s;
    
    k = k + 1;
    
    E_rb_all(:,:,:,k) = E_k;
end



%% Implementation of Eq 7a from Ulgubek Recursive Born
%  Using kth order born field computated previously
%  & forward propagating it to the image plane at z=0

% Multiply E with object function In spatial domain for every depth layer
% Eq(2) in Prof Tian's Compressive Bounds paper
%Ef=f.*E_k;

% Compute FFT of object-incident_field product
% (for filtering in fourier domain)
%F_Ef = F(f.*E_k);

% Multiply the object-incident_field FFT with phase3D in fourier domain
% & integrate along z-axis

cEsp = sum(F2(f.*E_k).*phase3D_H,3).*Pupil;

% Inverse fft to get scattered field at hologram plane
tmp = Ft2D(cEsp);
Holo_Efield_2D = tmp(1:Nx,1:Ny);

if isE2
    Holo_est = abs(Avg_amp+Holo_Efield_2D).^2;
else
    Holo_est = real(Holo_Efield_2D);
end

Perf_all = 0;
if born_order > 1
    born_orders = 1:born_order;
    for k = born_orders
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Implementation Performance metric evaluation
        %  Using kth order born field uk computated previously
        %  see how well it satisfies uk - G(uk.*f) = uin
        
        %1.Multiply kth order born field with the object
        Es=f.*E_rb_all(:,:,:,k);
        
        %2.Compute 2D FFTs for all slices of the previous product
        % Compute FFT of object-incident_field product
        % (for filtering in fourier domain)
        cEs = F(Es);
        cEsp =  cEs.*phase3D_G;
        temp= Ft3D(cEsp);
        E_G = temp(1:Nx,1:Ny,1:Nz);
        
        %4.Compute performance metric matrix
        perf_metric_mat = E_rb_all(:,:,:,k) - E_G - E;
        
        %5.Compute performance metric l2_norm_sq
        Perf_all(k) = norm(perf_metric_mat(:), 2)^2;
        
    end
end
end
