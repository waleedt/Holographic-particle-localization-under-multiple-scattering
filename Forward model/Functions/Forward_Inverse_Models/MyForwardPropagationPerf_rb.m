function perf_metric = MyForwardPropagationPerf_rb(f,E,Nx,Ny,Nz,...
    phase3D_G,born_order,F_nopad,F_pad,Ft,is_rm_alias)

if is_rm_alias
    cEs=zeros(3/2*Nx,3/2*Ny,Nz);
    F = F_pad;
else
    cEs=zeros(Nx,Ny,Nz);
    F = F_nopad;
end

% This matrix saves all born fields uk, uk-1, ... , u0
E_rb_all = zeros(Nx,Ny,Nz,born_order);
E_rb_all(:,:,:,1) = E;


%% Implementation of Eq 7b from Ulgubek Recursive Born
%  Recursive computation of kth order born field within the object
E_rb = E;
E_G = zeros(Nx,Ny,Nz);

k = 1;
while(k < born_order)
        
    % 1 Multiply k-1th order born field with the object
    Es=f.*E_rb;
    
    % Compute 2D FFTs for all slices of the previous product
    for i=1:Nz
        % Compute FFT of object-incident_field product
        % (for filtering in fourier domain)
        cEs(:,:,i)=F(Es(:,:,i));
    end
    
    % Now convolve with greens function of each depth separately
    % Propagate the eta.E_rb product, one by one, to all z-planes
    for j=1:Nz
        % Propagate eta.E_rb to the jth z-plane (first z-plane is
        % closest to z=0 image plane & j increases moving away from z=0
        %% Lei Tian comments:
        % the summation here is the discrete version of integration w.r.t.
        % to z, the integration size deltaZ is needed!
        %% Waleed comments:
        % deltaZ has been added in the definition of Phase3D_G
        cEsp=sum(cEs.*phase3D_G(:,:,:,j),3);
        
        temp= Ft(cEsp);
        E_G(:,:,j) = temp(1:Nx,1:Ny);
    end
    
    E_rb = E + E_G;
    
    k = k + 1;
    
    E_rb_all(:,:,:,k) = E_rb;
end


%% Implementation Performance metric evaluation
%  Using kth order born field uk computated previously 
%  see how well it satisfies uk - G(uk.*f) = uin

%1.Multiply kth order born field with the object
Es=f.*E_rb;

%2.Compute 2D FFTs for all slices of the previous product
for i=1:Nz
    % Compute FFT of object-incident_field product
    % (for filtering in fourier domain)
    cEs(:,:,i)=F(Es(:,:,i));
end

%3.Now convolve with greens function of each depth separately
% Propagate the eta.E_rb product, one by one, to all z-planes
for j=1:Nz
    % Propagate eta.E_rb to the jth z-plane (first z-plane is
    % closest to z=0 image plane & j increases moving away from z=0
    %% Lei Tian comments:
    % the summation here is the discrete version of integration w.r.t.
    % to z, the integration size deltaZ is needed!
    %% Waleed comments:
    % deltaZ has been added in the definition of Phase3D_G
    cEsp=sum(cEs.*phase3D_G(:,:,:,j),3);
    
    temp= Ft(cEsp);
    E_G(:,:,j) = temp(1:Nx,1:Ny);
end

%4.Compute performance metric matrix
perf_metric_mat = E_rb - E_G - E;

%5.Compute performance metric l2_norm_sq
perf_metric = norm(perf_metric_mat(:), 2)^2;
end
