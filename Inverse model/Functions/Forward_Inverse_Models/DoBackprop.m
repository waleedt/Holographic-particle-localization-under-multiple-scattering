% computes backpropagation (only in the 1st iteration)

% compute least sq solution by backpropagation
BPField = AT_init(Holo);
%BPField = reshape(MyV2C(f_least_sq),nx,ny,nz);

% Create folders for saving BPM results
bpm_dir = sprintf('%s/BPM',results_dir);
mkdir(bpm_dir);
mkdir([bpm_dir,'/BPM_abs']);
mkdir([bpm_dir,'/BPM_real']);
mkdir([bpm_dir,'/BPM_imag']);

% BPM Abs
fn =[bpm_dir,'/BPM_abs/','BPM_abs2_'];
bp_min= min(abs(BPField(:)));
bp_max = max(abs(BPField(:)));
for m =1 : nz
    %slice = mat2gray(abs(BPField(:,:,m)),[bp_min bp_max]);
    slice = uint8(((abs(BPField(:,:,m))- bp_max)/bp_min)*255);
    filename = [fn, num2str(m),'.png'];
    imwrite(slice,filename,'png')
end

% BPM Real
fn =[bpm_dir,'/BPM_real/','BPM_real2_'];
bp_min= min(real(BPField(:)));
bp_max = max(real(BPField(:)));
for m =1 : nz
    slice = mat2gray(real(BPField(:,:,m)),[bp_min bp_max]);
    filename = [fn, num2str(m),'.png'];
    imwrite(slice,filename,'png')
end

% BPM Imag
fn =[bpm_dir,'/BPM_imag/','BPM_imag2_'];
bp_min= min(imag(BPField(:)));
bp_max = max(imag(BPField(:)));

% Save BPM reconstruction as bmp images
for m =1 : nz
    slice = mat2gray(imag(BPField(:,:,m)),[bp_min bp_max]);
    filename = [fn, num2str(m),'.png'];
    imwrite(slice,filename,'png')
end
% no BPM for next iters
doBPM = 0;
% clean up memory
clear E0