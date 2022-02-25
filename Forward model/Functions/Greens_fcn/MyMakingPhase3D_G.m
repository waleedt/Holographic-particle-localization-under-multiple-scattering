function Phase3D_G = MyMakingPhase3D_G(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,F,is_rm_alias,isconv_G,isMemOpt)

%% Zero padding for ani aliasing
if is_rm_alias
    pad_len = 2;
    Nxs = Nx * 2;
    Nys = Ny * 2;
    Nx = Nx * pad_len;
    Ny = Ny * pad_len;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute via Weyl expansion
if(~isconv_G)
    if(~isMemOpt)
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
    else
        % code goes here
    end
    disp('phase3D_G_weyl done');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute via Conv approach
else
    k=(2*pi)/lambda;
    %% Spatial variables
    % x-y-z axes coordinates
    X_mat = kx * Nx * deltaX^2; % repeat X, 1 row Ny col
    Y_mat = ky * Ny * deltaY^2; % repeat Y, Nx row 1 col

    %% Generate Phase3D_G
    const = deltaX^2;%((k^2)/(4*pi))*
    
    if(~isMemOpt)
        Phase3D_G=zeros(Nx,Ny,Nz,Nz);
        for j=1:Nz
            for i=1:Nz
                Z_mat = repmat(abs(Z(j)-Z(i)),Nx,Ny);
                term = sqrt(X_mat.^2+Y_mat.^2+Z_mat.^2);
				term2D = sqrt(X_mat.^2+Y_mat.^2);
                if i~=j
                    Greens_spatial = exp(1i*k*term)./term;
                else
                    term(term==0)=0.1;
                    G_s_slice = exp(1i*k*term)./term;
                    G_s_slice(term2D<1) = (1/3)*(1/30);
                    Greens_spatial = G_s_slice;
                end
                Phase3D_G(:,:,i,j) = F(Greens_spatial)*const;
            end
        end
    else
        Z = [0:(Nz+1)-1].*deltaZ+offsetZ;
        Nzg = 2*(Nz+1)-1; %numslices in Phase3D_G
        Nzg_mid = (Nzg+1)/2; %idx of zero-propagation slice
        Phase3D_G=zeros(Nx,Ny,Nzg);
        for i=1:(Nz+1)
            Z_mat = repmat(abs(Z(1)-Z(i)),Nx,Ny);
            term = sqrt(X_mat.^2+Y_mat.^2+Z_mat.^2);
			term2D = sqrt(X_mat.^2+Y_mat.^2);
            if i==1
                term(term==0)=0.1; % avoid Nan
                Greens_spatial = exp(1i*k*term)./term;
                Greens_spatial(term2D<1) = (0.0116777 + 0.00330044i) * 1/30;
            else
                Greens_spatial = exp(1i*k*term)./term;
            end
            Phase3D_G(:,:,Nzg_mid-(i-1)) = F(Greens_spatial)*const;
            Phase3D_G(:,:,Nzg_mid+(i-1)) = F(Greens_spatial)*const;
            
        end
    end
    disp('phase3D_G_conv done');
end



