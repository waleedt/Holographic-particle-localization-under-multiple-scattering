% make output directories for results
born_order = forwardObj.bornorder;
f = xhat;
nz = forwardObj.sigSize(3);
iters_time = toc;

% compute normalized data fit
data_fit = (2*cost_fid)/(norm(Holo2D(:), 2)^2);

% Save f as mat file
fn =sprintf(results_dir_iter,iIter);
mkdir(fn);
filename =sprintf('%s/data.mat',fn);
if(isOracleSet)
    save(filename,'f','cost_all','cost_fid','cost_reg','iters_time',...
        'gradNorm','relativeChange_grad','z_est','snr_all','data_fit');
else
    save(filename,'f','cost_all','cost_fid','cost_reg','iters_time',...
        'gradNorm','relativeChange_grad','z_est','data_fit');
end

% Save hologram estimate
z_min= min(z_est(:));
z_max = max(z_est(:));
holo_img = mat2gray(z_est,[z_min,z_max]);
filename = sprintf('%s/holo_est.bmp',fn);
imwrite(holo_img,filename,'bmp')

% Save abs(f) as image files
ftmp = f;
% set max threshold
th1 = 0.7;
ftmp(ftmp>th1) = th1;
% scale to 0-255
ftmp = ftmp * (255/0.7);
% set binary threshold
th2 = 60;
ftmp(ftmp>=th2) = 255;
ftmp(ftmp<th2) = 0;
    
fig = figure(6);
for m =1 : nz
    imagesc(f(:,:,m)); axis image; colorbar;
    if m == 1
        c1 = caxis;
    else
        caxis(c1);
    end
    filename = sprintf('%s/TV_born%d_slice%d.png',fn,born_order,m);
    saveas(fig,filename);
end
filename = sprintf('%s/cost_snr.png',fn);
saveas(h,filename);
clear ftmp;