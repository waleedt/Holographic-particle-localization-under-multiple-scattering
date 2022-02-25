function grad=MyGradient(resid,f_temp,E_rb_all,Nx,Ny,Nz,phase3D_H,phase3D_G,...
    born_order,F_nopad,F_pad,Ft,is_rm_alias,Pupil,Avg_amp,SF,isE2,isMemOpt)

%%%
Nzg = 2*(Nz+1)-1; %numslices in Phase3D_G
Nzg_mid = (Nzg+1)/2; %idx of zero-propagation slice
%%%

if is_rm_alias
    cEs_sum = zeros(2*Nx,2*Ny,Nz);
    F = F_pad;
else
    cEs_sum = zeros(Nx,Ny,Nz);
    F = F_nopad;
end

% Compute gK = [Hh r] .* conj(uK)
%% Lei Tian comments:
% 1) the physical meaning of this step is to "back-propagate" the residual to 
%    the 3D object space, thus conj(.) is used.
% 2) the adjoint needs take into account the 2x factor in the 
%    hologram formation equation 2*real(E).
%    essentially, this is a quasi-Newton's method by taking into account
%    the non-unit svd values in the Forward model and try to fix it by
%    parts
if isE2
    resid_E2 = (Avg_amp+SF).*resid;
    cEsp=F(resid_E2).*Pupil;
else
    cEsp=F(resid).*Pupil;
end
%cEs=repmat(cEsp,[1 1 Nz]).*conj(phase3D_H);

% eta here denotes the virtual source due to the residual
if (phase3D_H == 1)
    tmp = Ft(repmat(cEsp,[1 1 Nz]).*conj(phase3D_G(:,:,Nzg_mid+1:end)));
else
    tmp = Ft(repmat(cEsp,[1 1 Nz]).*conj(phase3D_H));
end

eta = tmp(1:Nx,1:Ny,:);
clear tmp;

gK = conj(E_rb_all(:,:,:,born_order)).*eta;

% if born_order=3, we have K=2,g2, & E_rb corresponding to u2


if(born_order > 1)
    
    % Compute vK = [Hh r] .* conj(f)
    vK = f_temp.*eta;
    
    g_k_plus_1 = gK; %gK=g2
    v_k_plus_1 = vK;
    
    k = born_order;
    while (k > 1) %k=3 %k=2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Compute [Gh vk+1] .* conj(uk)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cEsp_3D = F(v_k_plus_1);
        
        for i = 1:Nz
            if(~isMemOpt)
                cEs_sum = cEs_sum + conj(squeeze(phase3D_G(:,:,i,:))).*repmat(cEsp_3D(:,:,i),[1 1 Nz]);
            else
                %Phase3D_G(:,:,:,3) = squeeze(Phase3D_G(:,:,3,:) checked
                cEs_sum = cEs_sum + conj(phase3D_G(:,:,Nzg_mid-(i-1):Nzg_mid-(i-1)+(Nz-1))).*repmat(cEsp_3D(:,:,i),[1 1 Nz]);
            end
        end
        
        tmp = Ft(cEs_sum);
        eta = tmp(1:Nx,1:Ny,:);
        clear tmp;
        
        %% Lei Tian commments:
        % need to consider the case when abs(E_rb_all) is not 1?
        %    Note: essentially, this is a quasi-Newton's method by taking into account
        %    the non-unit svd values in the Forward model and try to fix it by
        %    parts
        
        %born_order=3, we have k=3 for 1st iter, k-1 = 2, corresponding to u1  
        % E_rb corresponds to k-1=1 i.e. u0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Compute gk = gk+1 + [Gh vk+1] .* conj(uk)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        gk = g_k_plus_1 + eta.*conj(E_rb_all(:,:,:,k-1)); %g1 %g0

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Compute vk = [Gh vk+1] .* f
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %vk = f_temp.*eta; %v1 v0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Recursion of g & v variablezs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        g_k_plus_1 = gk;
        v_k_plus_1 = f_temp.*eta;    

        k = k - 1; %k = 2
    end
else
    gk = gK;
end

if isE2
    grad = 2*real(gk);
else
    grad = real(gk);
end

end