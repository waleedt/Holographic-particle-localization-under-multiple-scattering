function [Phase3D_H, Pupil]= MyMakingPhase3D_H_weyl(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,~,is_rm_alias,NA,ispupil)

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
kp=sqrt(kx.^2+ky.^2);
term=k.^2-kp.^2;
term(term<0)=0;

% z-coordinates of axial planes
Z = [0:Nz-1].*deltaZ+offsetZ;
Zcam = 0; % camera is at z=0


%% Generate Phase3D_H
Phase3D_H = zeros(Nx,Ny,Nz);

for i=1:Nz
    Phase3D_H(:,:,i) = (1i*pi*(k^2))*(exp(1i*2*pi*abs(Zcam-Z(i))*...
        sqrt(term))./(sqrt(term)))*deltaX;
end

%% Generate Pupil
if(ispupil)
    Pupil = kp <= (NA/lambda);   
else
    Pupil = ones(Nx,Ny);
end

end