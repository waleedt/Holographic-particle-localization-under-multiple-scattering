% save intermediate data for seagle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagate 3D field to image plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cEsp=sum(F_pad(f.*uk).*phase3D_H,3).*Pupil;
tmp = Ft(cEsp);
holo_efield_2D = tmp(1:Nx,1:Ny);

if isE2
    holo_seagle = abs(Avg_amp+holo_efield_2D).^2;
else
    holo_seagle = real(holo_efield_2D);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn =  sprintf('%s/iter_%d.mat',result_dir,k);
save(fn,'holo_seagle','holo_b1','holo_bm','uk','cost_vec',...
    'holo_efield_2D','f','E');