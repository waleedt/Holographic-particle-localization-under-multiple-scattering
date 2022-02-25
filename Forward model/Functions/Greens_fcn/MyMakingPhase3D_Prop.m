function [Phase3D_Prop, Pupil_prop, k2_bpm] = MyMakingPhase3D_Prop(Nx,Ny,Nz,...
    lambda,deltaX,deltaY,deltaZ,offsetZ,NA,ispupil)

% Spatial frequency variables
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

% Generate Phase3D_Prop
Phase3D_Prop=zeros(Nx,Ny,Nz);
for i=1:Nz
    %% Lei Tian Comments:
    % Angular spectrum representation of forward propagator, which is different
    % from the Green's function angular spectrum expansion (also known as
    Phase3D_Prop(:,:,i)=exp(1j*2*pi*(Zcam-Z(i))*sqrt(term));
end

% Generate Pupil
if(ispupil)
    Pupil_prop = kp <= (NA/lambda);
    Pupil_prop = fftshift(Pupil_prop);
else
    Pupil_prop = 1;
end

k2_bpm = sqrt(term);

end