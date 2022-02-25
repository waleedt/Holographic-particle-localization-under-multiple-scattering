function Phase3D_G = MyMakingPhase3D_G_weyl(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,~,is_rm_alias,~)

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