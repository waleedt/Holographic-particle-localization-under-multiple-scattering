function numgrad=NumGradient(holo,f,eps,E,born_order,phase3D_H,phase3D_G,numiters)

[Nx,Ny,Nz] = size(f);
N = length(f(:));
numgrad = zeros(Nx,Ny,Nz);

if numiters ~= 0
    N = numiters;
end

for iii = 1:N
    % convert linear index to 3D index
    siz = size(f);
    [ii,jj,kk] = ind2sub(siz,iii);
    
    %% compute real part of gradient first
    % compute J+ ans J-
    ei = zeros(Nx,Ny,Nz);
    ei(ii,jj,kk) = 1;
    
    f_plus = f + eps * ei;
    f_minus = f - eps * ei;
    
    % J+
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Recursive computation of kth order born field within the object
    E_k = E;
    E_s  = E;
    k = 1;
    while(k < born_order)
        % 1 Multiply k-1th order born scatteredfield with the object
        Ef=E_s.*f_plus;
        % Compute 2D FFTs for all slices of the previous product
        F_Ef = zeros(2*Nx,2*Ny,Nz);
        for i=1:Nz
            F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),2*size(Ef(:,:,i),2)));
        end
        
        % Now convolve with greens function of each depth separately
        for j=1:Nz
            cEsp=sum(F_Ef.*phase3D_G(:,:,:,j),3);
            tmp = ifft2(ifftshift(cEsp));
            E_s(:,:,j) = tmp(1:Nx,1:Ny);
        end
        E_k = E_k + E_s;
        k = k + 1;
    end
    %  Using kth order born field computated previously to propagating it to the image plane
    Ef=f_plus.*E_k;
    for i=1:Nz
        F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),2*size(Ef(:,:,i),2)));
    end
    % Multiply the object-incident_field FFT with phase3D in fourier domain & integrate along z-axis
    cEsp=sum(F_Ef.*phase3D_H,3);
    % Inverse fft to get scattered field at hologram plane
    tmp = ifft2(ifftshift(cEsp));
    Holo_Efield_2D = tmp(1:Nx,1:Ny);
    % Hologram is Re(scattered_field) at hologram plane as per derivation
    holo_est = real(Holo_Efield_2D);
    resid = holo_est - holo;
    J_plus = 0.5*norm(resid(:),2)^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % J-
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Recursive computation of kth order born field within the object
    E_k = E;
    E_s  = E;
    k = 1;
    while(k < born_order)
        % 1 Multiply k-1th order born scatteredfield with the object
        Ef=E_s.*f_minus;
        % Compute 2D FFTs for all slices of the previous product
        for i=1:Nz
            F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),2*size(Ef(:,:,i),2)));
        end
        
        % Now convolve with greens function of each depth separately
        for j=1:Nz
            cEsp=sum(F_Ef.*phase3D_G(:,:,:,j),3);
            tmp = ifft2(ifftshift(cEsp));
            E_s(:,:,j) = tmp(1:Nx,1:Ny);
        end
        E_k = E_k + E_s;
        k = k + 1;
    end
    %  Using kth order born field computated previously to propagating it to the image plane
    Ef=f_minus.*E_k;
    for i=1:Nz
        F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),2*size(Ef(:,:,i),2)));
    end
    % Multiply the object-incident_field FFT with phase3D in fourier domain & integrate along z-axis
    cEsp=sum(F_Ef.*phase3D_H,3);
    % Inverse fft to get scattered field at hologram plane
    tmp = ifft2(ifftshift(cEsp));
    Holo_Efield_2D = tmp(1:Nx,1:Ny);
    % Hologram is Re(scattered_field) at hologram plane as per derivation
    holo_est = real(Holo_Efield_2D);
    resid = holo_est - holo;
    J_minus = 0.5*norm(resid(:),2)^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % copmtue gradient value
    grad_re = (J_plus - J_minus) / (2*eps);
    
    %% compute imaginary part of gradient
    % compute J+ ans J-
    f_plus = f + eps * ei * 1i;
    f_minus = f - eps * ei* 1i;
    
    % J+
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Recursive computation of kth order born field within the object
    E_k = E;
    E_s  = E;
    k = 1;
    while(k < born_order)
        % 1 Multiply k-1th order born scatteredfield with the object
        Ef=E_s.*f_plus;
        % Compute 2D FFTs for all slices of the previous product
        F_Ef = zeros(2*Nx,2*Ny,Nz);
        for i=1:Nz
            F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),2*size(Ef(:,:,i),2)));
        end
        
        % Now convolve with greens function of each depth separately
        for j=1:Nz
            cEsp=sum(F_Ef.*phase3D_G(:,:,:,j),3);
            tmp = ifft2(ifftshift(cEsp));
            E_s(:,:,j) = tmp(1:Nx,1:Ny);
        end
        E_k = E_k + E_s;
        k = k + 1;
    end
    %  Using kth order born field computated previously to propagating it to the image plane
    Ef=f_plus.*E_k;
    for i=1:Nz
        F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),2*size(Ef(:,:,i),2)));
    end
    % Multiply the object-incident_field FFT with phase3D in fourier domain & integrate along z-axis
    cEsp=sum(F_Ef.*phase3D_H,3);
    % Inverse fft to get scattered field at hologram plane
    tmp = ifft2(ifftshift(cEsp));
    Holo_Efield_2D = tmp(1:Nx,1:Ny);
    % Hologram is Re(scattered_field) at hologram plane as per derivation
    holo_est = real(Holo_Efield_2D);
    resid = holo_est - holo;
    J_plus = 0.5*norm(resid(:),2)^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % J-
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Recursive computation of kth order born field within the object
    E_k = E;
    E_s  = E;
    k = 1;
    while(k < born_order)
        % 1 Multiply k-1th order born scatteredfield with the object
        Ef=E_s.*f_minus;
        % Compute 2D FFTs for all slices of the previous product
        for i=1:Nz
            F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),2*size(Ef(:,:,i),2)));
        end
        
        % Now convolve with greens function of each depth separately
        for j=1:Nz
            cEsp=sum(F_Ef.*phase3D_G(:,:,:,j),3);
            tmp = ifft2(ifftshift(cEsp));
            E_s(:,:,j) = tmp(1:Nx,1:Ny);
        end
        E_k = E_k + E_s;
        k = k + 1;
    end
    %  Using kth order born field computated previously to propagating it to the image plane
    Ef=f_minus.*E_k;
    for i=1:Nz
        F_Ef(:,:,i)=fftshift(fft2(Ef(:,:,i),2*size(Ef(:,:,i),1),2*size(Ef(:,:,i),2)));
    end
    % Multiply the object-incident_field FFT with phase3D in fourier domain & integrate along z-axis
    cEsp=sum(F_Ef.*phase3D_H,3);
    % Inverse fft to get scattered field at hologram plane
    tmp = ifft2(ifftshift(cEsp));
    Holo_Efield_2D = tmp(1:Nx,1:Ny);
    % Hologram is Re(scattered_field) at hologram plane as per derivation
    holo_est = real(Holo_Efield_2D);
    resid = holo_est - holo;
    J_minus = 0.5*norm(resid(:),2)^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % copmtue gradient value
    grad_im = (J_plus - J_minus) / (2*eps*1i);
    
    %% compute total gradient
    numgrad(iii) = grad_re + 1i*grad_im;
    
    fprintf('iter:%d\n',iii);
end

end