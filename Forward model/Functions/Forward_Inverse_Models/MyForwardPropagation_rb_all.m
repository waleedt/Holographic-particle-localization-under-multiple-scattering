% compute holograms corresponding to all the born orders listed in the
% born_orders array
function [Holo_px_sampled,Holo_not_px_sampled,Holo_Efield_2D,E_rb_all,perf_metric] = ...
    MyForwardPropagation_rb_all(f,E,Nx,Ny,Nz,phase3D_H,phase3D_G,...
    born_orders,F_nopad,F_pad,Ft,is_rm_alias,is_perf_metric,deltaZ,...
    super_res_factor)

%
if is_rm_alias
    F_Ef=zeros(2*Nx,2*Ny,Nz);
    F = F_pad;
else
    F_Ef=zeros(Nx,Ny,Nz);
    F = F_nopad;
end

Holo_not_px_sampled = zeros(Nx,Ny,length(born_orders));
Holo_px_sampled = zeros(Nx/super_res_factor,Ny/super_res_factor,length(born_orders));
Holo_Efield_2D = zeros(Nx,Ny,length(born_orders));
perf_metric = zeros(2,length(born_orders));

% get the max born order in the born_orders array
max_born_order = max(born_orders);

% This matrix saves all born fields uk, uk-1, ... , u0
E_rb_all = zeros(Nx,Ny,Nz,max_born_order);
E_rb_all(:,:,:,1) = E;

%% Implementation of Eq 7b from Ulgubek Recursive Born
%  Recursive computation of kth order born field within the object
E_k = E;
E_s  = E;

k = 1;
while(k < max_born_order)
        
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


%% Now that we have all multiple scattered born fields in the object,
% in order to compute the hologram corresponding to a certain born-order,
% we will simply pointwise multiply that born field with the object and
% bpropagate it to the hologram plane

% Implementation of Eq 7a from Ulgubek Recursive Born
%  Using kth order born field computated previously 
%  & forward propagating it to the image plane at z=0

for k = 1:length(born_orders)
    % Multiply E with object function In spatial domain for every depth layer
    % Eq(2) in Prof Tian's Compressive Bounds paper
    Ef=f.*E_rb_all(:,:,:,born_orders(k));

    for i=1:Nz
        % Compute FFT of object-incident_field product
        % (for filtering in fourier domain)
        %F_Ef(:,:,i)=F(Ef(:,:,i));
        F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),...
                2*size(Ef(:,:,i),2))); % above line hardcoded
    end

    % Multiply the object-incident_field FFT with phase3D in fourier domain
    % & integrate along z-axis
    cEsp=sum(F_Ef.*phase3D_H,3);

    % Inverse fft to get scattered field at hologram plane
    tmp = Ft(cEsp);
    Holo_Efield_2D(:,:,k) = tmp(1:Nx,1:Ny);

    % Hologram is 2*scattered_field at hologram plane
    Holo_not_px_sampled(:,:,k) = 1 + 2*real(Holo_Efield_2D(:,:,k));

    Holo_px_sampled(:,:,k) = pixel_sampling(Holo_not_px_sampled(:,:,k),super_res_factor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Metric Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    perf_metric(1,k) = norm(perf_metric_mat(:), 2);
    % compute norm of 3D scattered field
    E_s_k = zeros(Nx,Ny,Nz);
    if born_orders(k) > 1
        E_s_k = E_rb_all(:,:,:,born_orders(k)) - E_rb_all(:,:,:,born_orders(k)-1);
    end
    perf_metric(2,k) = norm(E_s_k(:), 2);
else
    perf_metric(1,k) = NaN;
end
end

end
