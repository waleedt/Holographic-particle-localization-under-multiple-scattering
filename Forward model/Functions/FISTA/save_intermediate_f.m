% make output directories for results
born_order = forwardObj.bornorder;
f = xhat;
nz = forwardObj.sigSize(3);
iters_time = toc;

f_outdir = sprintf('%s/born%d/iters_%d',results_dir,born_order,iIter);
mkdir(f_outdir);
mkdir([f_outdir,'/f_img']);

% compute normalized data fit
data_fit = (2*cost_fid)/(norm(Holo2D(:), 2)^2);

%%% f-reconstruction result - saving

% Save f as mat file
fn =[f_outdir,'/f_born',num2str(born_order),'.mat'];
if(isOracleSet)
    save(fn,'f','cost_all','cost_fid','cost_reg','iters_time',...
        'gradNorm','relativeChange_grad','z_est','snr_all','data_fit');
else
    save(fn,'f','cost_all','cost_fid','cost_reg','iters_time',...
        'gradNorm','relativeChange_grad','z_est','data_fit');
end

% Save hologram estimate
z_min= min(z_est(:));
z_max = max(z_est(:));
holo_img = mat2gray(z_est,[z_min,z_max]);
fn_h =[f_outdir,'/holo_est.bmp'];
imwrite(holo_img,fn_h,'bmp')

% Save abs(f) as image files
fn =sprintf(results_dir_iter,iIter);
mkdir(fn);
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
    
for m =1 : nz
    slice = mat2gray(ftmp(:,:,m));
    filename = [fn,'/TV_born',num2str(born_order),'_slice_', num2str(m),'.png'];
    imwrite(slice,filename,'png');
end
clear ftmp;