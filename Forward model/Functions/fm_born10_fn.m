% n_obj = 1.52 or 1.34 for high or low contrast
function fm_born10_fn(n_obj,resfolder)
add_paths; func_handles;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Parameters - Primary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stype = 3; %for 3D Born

% Rg to par_per_slice map
obj_slices = containers.Map({1,2,3}, {1,2,3});
obj_stype = containers.Map({1,2,3}, {'disc_const','disc_int','sphere'});
obj_conc = containers.Map({1,2,4,8,16}, {1,2,4,8,16});

% load data
filename = sprintf('%s/f_obj_1.mat',resfolder);
load(filename);
f = f_obj;
[nx,ny,rec_num_slices] = size(f);

%%%%% Forward model parameter settings (fixed)
%%%% BPM
% remove circular convolution aliasing during propagation
remalias_bpm = 1;
% is there paraxial approximation in the model
isparaxial_bpm = 0;
%%% RBorn 
% E2_incl in hologram
isE2 = 1;
% super resolution in object domain w.r.t pixel size
super_res_factor = 1;
% remove aliasing in convolutions
is_rm_alias = 1;
% copmute performance metric
is_perf_metric = 1;
% limit NA by pupil fcn
ispupil = 0;
% method of computing Phase3D_G [0-Weyl, 1-Conv]
isconv_G = 1;
% method of computing Phase3D_H [0-Weyl, 1-Conv]
isconv_H = 1;
% optimize memory usage?
isMemOpt = 1;

% Incident field amplitude and intensity
Avg_int = 1;
Avg_amp = sqrt(1);

% initialize simulation parameters
sim_params_3_45_20x;
f = f * (n_medium^2 - n_obj^2);

%%% Initialize incident field and propagation operators
% Define incident field magnitude
E0  = ones(nx,ny).*Avg_amp;
% Phase3D_H defines propagation from 3D object to 2D hologram
[Phase3D_H, Pupil, k_bpm]= MyMakingPhase3D_H(nx,ny,nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,FG,is_rm_alias,isconv_H,ispupil,NA,...
    remalias_bpm,isparaxial_bpm);
% Phase3D_G defines propagation from 3D object to within the object
Phase3D_G = MyMakingPhase3D_G_3D(nx,ny,nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,FG3D,is_rm_alias,isconv_G,isMemOpt);
% Phase3D_Prop defines unhindered propagation of E0 within the object
[Phase3D_Prop, Pupil_prop, k2_bpm] = MyMakingPhase3D_Prop(nx,ny,nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,NA,ispupil);
% Back-propagating E0 in free space from z=0 to within obj cube
E = MyFieldsBackPropagation(E0,nx,ny,nz,Phase3D_Prop,F,Ft,Pupil_prop);
BPM_E0slice = E(:,:,end);

% create results directory
result_dir = sprintf('%s/n%1.2f/',resfolder,n_obj);
mkdir(result_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Born
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_rb = @(f,born_order) MyForwardPropagation_rb3D(f,E,nx,ny,nz,Phase3D_H,...
    Phase3D_G,born_order,F3D,Fpad3D,F,Fpad,Ft3D,Ft,is_rm_alias,...
    Pupil,Avg_amp,isE2,isMemOpt);

[holo_rb3D,~,~,perf] = A_rb(f,15);

figure;imagesc(real(holo_rb3D)); axis image; colorbar; colormap(hot); title('rBorn');

figure;
plot(perf,'-*','linewidth',2);
xlabel('Born order');
ylabel('Perf metric');
title('Convergence of Born');

fn = sprintf('%s/holo_rb3D.mat',result_dir);
save(fn,'holo_rb3D');

end
