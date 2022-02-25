function [Holo_est,perf_metric]=...
    MyForwardPropE2_rb(f,E,Nx,Ny,Nz,phase3D_H,phase3D_G,...
    born_order,F_nopad,F_pad,Ft,is_rm_alias,Avg_amp,Pupil)

%
if is_rm_alias
    F_Ef=zeros(2*Nx,2*Ny,Nz);
    F = F_pad;
else
    F_Ef=zeros(Nx,Ny,Nz);
    F = F_nopad;
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
    Ef=E_s.*f;
    
    % Compute 2D FFTs for all slices of the previous product
    for i=1:Nz
        % Compute FFT of object-incident_field product
        % (for filtering in fourier domain)
        
        %F_Ef(:,:,i)=F(Ef(:,:,i));
        F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),...
            2*size(Ef(:,:,i),2))); % above line hardcoded
    end
    
    % Now convolve with greens function of each depth separately
    % Propagate the eta.E_k product, one by one, to all z-planes
    for j=1:Nz
        % Propagate eta.E_k to the jth z-plane (first z-plane is
        % closest to z=0 image plane & j increases moving away from z=0
        %% Lei Tian comments:
        % the summation here is the discrete version of integration w.r.t.
        % to z, the integration size deltaZ is needed!
        %% Waleed comments:
        % deltaZ has been added in the definition of Phase3D_G
        %cEsG = F_Ef;
        %cEsG(:,:,j) = F(E_k(:,:,j)); %added
        cEsp=sum(F_Ef.*phase3D_G(:,:,:,j),3);
        
        %tmp = Ft(cEsp);
        tmp = ifft2(ifftshift(cEsp)); % above line hardcoded
        
        E_s(:,:,j) = tmp(1:Nx,1:Ny);
    end
    
    E_k = E_k + E_s;
    
    k = k + 1;
    
    E_rb_all(:,:,:,k) = E_k;
end


%% Implementation of Eq 7a from Ulgubek Recursive Born
%  Using kth order born field computated previously 
%  & forward propagating it to the image plane at z=0

% Multiply E with object function In spatial domain for every depth layer
% Eq(2) in Prof Tian's Compressive Bounds paper
Ef=f.*E_k;

for i=1:Nz
    % Compute FFT of object-incident_field product
    % (for filtering in fourier domain)
    %F_Ef(:,:,i)=F(Ef(:,:,i));
    F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),...
            2*size(Ef(:,:,i),2))); % above line hardcoded
end

% Multiply the object-incident_field FFT with phase3D in fourier domain
% & integrate along z-axis
cEsp=sum(F_Ef.*phase3D_H,3).*Pupil;

% Inverse fft to get scattered field at hologram plane
tmp = Ft(cEsp);
Holo_Efield_2D = tmp(1:Nx,1:Ny);

Holo_est = abs(Avg_amp+Holo_Efield_2D).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Metric Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_perf_metric = 1;
if is_perf_metric
    %steps 1 & 2 are in lines 67 & 69-73 respectively
    %step 3: Now convolve with greens function of each depth separately
    % Propagate the eta.E_k product, one by one, to all z-planes
    for j=1:Nz
        % Propagate E_k.*f to the jth z-plane (first z-plane is
        % closest to z=0 image plane & j increases moving away from z=0
        cEsp=sum(F_Ef.*phase3D_G(:,:,:,j),3);
        
        temp= Ft(cEsp);
        E_s(:,:,j) = temp(1:Nx,1:Ny);
    end
    
    %step 4: Compute performance metric matrix
    perf_metric_mat = E_k - E_s - E;
    
    %step 5: Compute performance metric l2_norm
    perf_metric = norm(perf_metric_mat(:), 2);
else
    perf_metric = NaN;
end