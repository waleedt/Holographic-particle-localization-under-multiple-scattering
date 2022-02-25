function [grad, b] = Grad_seagle(u,uin,f,phase3D_G,Nx,Ny,Nz,F_pad,Ft)

%%%
Nzg = 2*(Nz+1)-1; %numslices in Phase3D_G
Nzg_mid = (Nzg+1)/2; %idx of zero-propagation slice
%%%

cEs_sum = zeros(2*Nx,2*Ny,Nz);
G_diagF_u = zeros(Nx,Ny,Nz);

F = F_pad;


%% compute b
diagF_u_fft = F(f.*u);

% filtering with G
for j=1:Nz
    cEsp=sum(diagF_u_fft.*phase3D_G(:,:,Nzg_mid-(j-1):Nzg_mid-(j-1)+(Nz-1)),3);
    
    temp= Ft(cEsp);
    G_diagF_u(:,:,j) = temp(1:Nx,1:Ny);
end

b = u - G_diagF_u - uin;

%% compute gradient

% compute Gh_b

cEsp_3D = F(b);

for i = 1:Nz
    cEs_sum = cEs_sum + conj(phase3D_G(:,:,Nzg_mid-(i-1):Nzg_mid-(i-1)+(Nz-1))).*repmat(cEsp_3D(:,:,i),[1 1 Nz]);
end

tmp = Ft(cEs_sum);
Gh_b = tmp(1:Nx,1:Ny,:);
clear tmp;

% gradient
grad = b -  f .* Gh_b;

end