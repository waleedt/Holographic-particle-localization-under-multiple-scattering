function f_lest_sq=MyAdjointPropagation(Holo2D,E,Nx,Ny,Nz,phase3D,...
    F,Ft,sv_sq,Pupil)

cEsp=F(Holo2D).*conj(Pupil);

cEs=repmat(cEsp,[1 1 Nz]).*conj(phase3D);

invf_cEs=zeros(Nx,Ny,Nz);
for i=1:Nz
    tmp = Ft(cEs(:,:,i));
    invf_cEs(:,:,i)=tmp(1:Nx,1:Ny);
end

%% Lei Tian commments:
% need to consider the case when abs(E) is not 1
f_lest_sq=conj(E).*invf_cEs;

%% Waleed comments 1/17/17
% Ah(holo) = sv^2 * eta
% eta = Ah(holo) * 1/sv^2
f_lest_sq = f_lest_sq * 1/sv_sq;
end
