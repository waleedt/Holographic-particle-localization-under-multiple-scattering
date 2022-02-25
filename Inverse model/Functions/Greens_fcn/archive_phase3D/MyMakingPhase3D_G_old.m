function Phase3D_G = MyMakingPhase3D_G_old(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,F,is_rm_alias,isconv_G)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute via Weyl expansion
if(~isconv_G) 
    %% Zero padding for ani aliasing
    if is_rm_alias
        Nx = Nx * 2;
        Ny = Ny * 2;
    end
    
    %% Spatial frequency variables
    k=1/lambda;
    KX=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));
    KY=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));
    kx=repmat(KX,1,Ny); % repeat X, 1 row Ny col
    ky=repmat(KY,Nx,1); % repeat Y, Nx row 1 col
    kp=kx.^2+ky.^2;
    term=k.^2-kp;
    term(term<0)=0;
    
    % z-coordinates of axial plane
    Z = [0:Nz-1].*deltaZ+offsetZ;
    
    %% Generate Phase3D_G
    Phase3D_G = zeros(Nx,Ny,Nz,Nz);
    for j=1:Nz % separate Phase3D for each z
        for i=1:Nz
            if i~=j
                % propagation from slice at z(i) to sclice at z(j)
                Phase3D_G(:,:,i,j) = (1j*pi*(k^2)*deltaX)*exp(1j*2*pi*abs(Z(j)-Z(i))*sqrt(term))./(sqrt(term));
            else
                Phase3D_G(:,:,i,j)=zeros(Nx,Ny);
            end
            
        end
    end
    disp('phase3D_G_weyl done');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else % Compute via Conv approach
    %% Zero padding for ani aliasing
    if is_rm_alias
        pad_len = 2;
        Nxs = Nx * 2;
        Nys = Ny * 2;
        Nx = Nx * pad_len;
        Ny = Ny * pad_len;
    end
    
    %% Spatial variables
    % x-y-z axes coordinates
    k=(2*pi)/lambda;
    X = [ceil(-Nx/2):1:ceil(Nx/2-1)]'*deltaX;
    Y = [ceil(-Ny/2):1:ceil(Ny/2-1)]*deltaY;
    Z = [0:Nz-1].*deltaZ+offsetZ;
    X_mat=repmat(X,1,Ny); % repeat X, 1 row Ny col
    Y_mat=repmat(Y,Nx,1); % repeat Y, Nx row 1 col
    % calculating phase shift in x and y
    kx = 0:Nx-1;
    ky = 0:Ny-1;
    phase_shiftx = exp(+1i*(2*pi*kx/(Nx))*(Nx/2));
    phase_shifty = exp(-1i*(2*pi*ky/(Ny))*(Ny/2));
    Phase_x = repmat(phase_shiftx,Nx,1); % FFT compensated for phase shift in  x
    Phase_y = repmat(phase_shifty',1,Ny); % FFT compensated for phase shift in  y
    
    %% Generate Phase3D_G
    Phase3D_G=zeros(Nx,Ny,Nz,Nz);
    for j=1:Nz
        for i=1:Nz
            if i~=j
                Z_mat = repmat(abs(Z(j)-Z(i)),Nx,Ny);
                term = sqrt(X_mat.^2+Y_mat.^2+Z_mat.^2);
                const = ((k^2)/(4*pi))*deltaX^3;
                Phase3D_G(:,:,i,j) = (F(exp(1i*k*term)./term).*Phase_x.*Phase_y)*const;
            else
                Phase3D_G(:,:,i,j)=zeros(Nxs,Nys);
            end
        end
    end
    disp('phase3D_G_conv done');
end

