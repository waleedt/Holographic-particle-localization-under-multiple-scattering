function [Phase3D_H, Pupil]= MyMakingPhase3D_H(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,F,is_rm_alias,isconv_H,ispupil,NA)

%% Zero padding for ani aliasing
if is_rm_alias
    pad_len = 2;
    Nx = Nx * pad_len;
    Ny = Ny * pad_len;
end

%% Spatial frequency variables
k=1/lambda;
KX=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));
KY=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));
kx=repmat(KX,1,Ny); % repeat X, 1 row Ny col
ky=repmat(KY,Nx,1); % repeat Y, Nx row 1 col
kp=sqrt(kx.^2+ky.^2);
term=k.^2-kp.^2;
term(term<0)=0;

% z-coordinates of axial planes
Z = [0:Nz-1].*deltaZ+offsetZ;
Zcam = 0; % camera is at z=0

%% Generate Pupil
    if(ispupil)
        Pupil = double(kp <= (NA/lambda));
        Pupil = fftshift(Pupil);
    else
        Pupil = ones(Nx,Ny);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute via Weyl expansion
if(~isconv_H)
    %% Generate Phase3D_H
    Phase3D_H = zeros(Nx,Ny,Nz);
    
    for i=1:Nz
        Phase3D_H(:,:,i) = (1i*pi*(k^2))*(exp(1i*2*pi*abs(Zcam-Z(i))*...
            sqrt(term))./(sqrt(term)))*deltaX;
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute via Conv approach
else 
    k = 2*pi/lambda;
    %% Spatial variables
    % x-y axes coordinates
    X_mat = kx * Nx * deltaX^2; % repeat X, 1 row Ny col
    Y_mat = ky * Ny * deltaY^2; % repeat Y, Nx row 1 col
    
    %% Generate Phase3D_H
    Phase3D_H=zeros(Nx,Ny,Nz);
    Greens_spatial = zeros(Nx,Ny,Nz);
    dist = zeros(Nx,Ny,Nz);
    debg = 0; % debugging
    for i=1:Nz
        Z_mat = repmat(Z(i),Nx,Ny);
        term = sqrt(X_mat.^2+Y_mat.^2+Z_mat.^2);
        if(debg)
            dist(:,:,i) = sqrt(X_mat.^2+Y_mat.^2+Z_mat.^2);
            Greens_spatial(:,:,i) = exp(1i*k*term)./term;
        end
        % shifting because matlab starts counting from 0, whereas, we want
        % the counting to start form -ve. Multiplying with const term k^2/4pi
        % for object scattering density and dV for the area under unit impulse
        const = ((k^2)/(4*pi))*deltaX^3;
        Phase3D_H(:,:,i) = F(exp(1i*k*term)./term)*const;
    end
end

end