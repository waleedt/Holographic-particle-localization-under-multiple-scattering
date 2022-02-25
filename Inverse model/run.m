
clc; clear all; close all;
add_paths; func_handles;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Parameters - Primary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% born order for reconstruction
rec_order = 1;
% object ref index
n_obj = 1.43;
% reg parameter
tau_ = 0.07;%[0.05 0.07 0.1 0.12 0.15 0.2];

% px spacing of object
for px_spacing = 8:8
    
    fprintf('************** px_spacing_%d ****************',px_spacing);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set Parameters - Secondary
    % numslices for reconstruction
    rec_num_slices = 3;
    % lateral size of reconstruction
    holo_size = 256;
    
    % E2_incl?
    isE2 = 0;
    % do BPM?
    doBPM = 1;
    % save_every x iters
    save_every_f = 100;
    
    % LOAD DATA
    fn = sprintf('Data/px_spacing_%d/n%1.2f/holo_rb3D.mat',px_spacing,n_obj);
    load(fn);
    Holo = (holo_rb3D - 1) / 2;
    fn = sprintf('Data/px_spacing_%d/f_obj_1.mat',px_spacing);
    load(fn);
    f_gt = sum(f_obj(:,:,1:8),3);
    if (rec_num_slices == 3)
        f_gt(:,:,2) = zeros(256,256);
        f_gt(:,:,3) = f_gt(:,:,1);
    end
    
    super_res_factor = 1; % super resolution in object domain w.r.t pixel size
    is_rm_alias = 1;      % remove aliasing in convolutions? currenlty the code
    %only works when this is set to 1. Hardcoding is in the parfor loops
    is_perf_metric = 0;
    % limit NA by pupil fcn - flag
    ispupil = 0;
    % method of computing Phase3D_G [0-Weyl, 1-Conv]
    isconv_G = 1;
    % method of computing Phase3D_H [0-Weyl, 1-Conv]
    isconv_H = 1;
    % optimize memory usage?
    isMemOpt = 1;
    
    %%% Main Computation Loop
    for tau = tau_
        for born_order = rec_order
            
            tic
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Hologram preprocessing
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % initialize simulation parameters
            nx = holo_size;
            ny = holo_size;
            sim_params_3_45_20x;
            % Simulate object
            %f_gt = sim_4circ_three_6px(nx,ny,nz,n_obj,n_medium,1);
            f_gt = f_gt * (n_medium^2 - n_obj^2);
            % incident field amp and int
            Avg_int = 1;
            Avg_amp = sqrt(Avg_int);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % initialize incident field and propagation operators
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % define incident field magnitude
            % Lei Tian edits:
            % mean2(Holo) is the estimated DC intensity of input field
            E0  = ones(nx,ny).*Avg_amp;
            
            % Phase3D_H defines propagation from 3D object to 2D hologram
            [Phase3D_H, Pupil]= MyMakingPhase3D_H(nx,ny,nz,lambda,...
                deltaX,deltaY,deltaZ,offsetZ,FG,is_rm_alias,isconv_H,ispupil,NA);
            
            % Phase3D_G defines propagation from 3D object to within the object
            % (repetitive convolution with this yeilds higher order fields)
            Phase3D_G = MyMakingPhase3D_G(nx,ny,nz,lambda,...
                deltaX,deltaY,deltaZ,offsetZ,FG,is_rm_alias,isconv_G,isMemOpt);
            % Phase3D_Prop defines unhindered propagation of E0 within the object
            % Lei Tian comments:
            % field propagatino needs to be done based on Angular spectrum,
            % which is slightly different from the Green's function.
            [Phase3D_Prop, Pupil_prop] = MyMakingPhase3D_Prop(nx,ny,nz,lambda,...
                deltaX,deltaY,deltaZ,offsetZ,NA,ispupil);
            
            % back-propagating E0 in free space from z=0 to within obj cube
            E = MyFieldsBackPropagation(E0,nx,ny,nz,Phase3D_Prop,F,Ft,Pupil_prop);
            clear Phase3D_Prop;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Function handles for forward (A) and inverse (AT) oprators
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % forward operator
            A = @(f) MyForwardPropagation_rb(f,E,nx,ny,nz,Phase3D_H,...
                Phase3D_G,born_order,F,Fpad,Ft,is_rm_alias,...
                Pupil,Avg_amp,isE2,isMemOpt);
            
            Asim = @(f) MyForwardPropagation_rb_all(f,E,nx,ny,nz,Phase3D_H,...
                Phase3D_G,born_orders,F,Fpad,Ft,is_rm_alias,...
                Pupil,Avg_amp,isE2,isMemOpt);
            
            % compute singular values of A
            % singular values of A are sqrt(eig(A'A))
            % where, eig(A'A) = 2*(Avg_amp*nz*pi^2*deltaZ^2)/lambda^2
            sv_sq = (Avg_int*nz*pi^2*deltaZ^2)/lambda^2;
            
            
            % adjoint/backward operator
            Grad = @(resid,f,E_rb_all,SF) MyGradient(resid,f,E_rb_all,nx,ny,nz,...
                Phase3D_H,Phase3D_G,born_order,F,Fpad,Ft,is_rm_alias,Pupil,...
                Avg_amp,SF,isE2,isMemOpt);
            
            % AT_init does not depend on f, AT does. Used for initial guess of f
            AT_init = @(g) MyAdjointPropagation(g,E,nx,ny,nz,Phase3D_H,Fpad,...
                Ft,sv_sq,Pupil);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % forward model simulation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %         [Holos_all,~,~,perf_all] = Asim(f_gt);
            %         Holo = Holos_all(:,:,end);
            %         figure(2)
            %         plot(perf_all,'-*','linewidth',2);
            %         grid on;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % backpropagation reconstruction (least sq solution)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            results_dir = sprintf('Results_iter_%1.2f/px_spacing_%d/',n_obj,px_spacing);
            results_dir_iter = sprintf('Results_iter_%1.2f/px_spacing_%d/tau_%1.0e/iter%%d/born%d',n_obj,px_spacing,tau,born_order);
            mkdir(results_dir)
            
            if doBPM
                DoBackprop;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Run FISTA to get final result. Ulugbek's implementation 'FPPA'.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% Generate forward model
            fwd_inv_obj = FwdInvModel(A,Grad,AT_init,born_order,[nx,ny,nz]);
            
            %%% Run FISTA
            [f, ~] = fistaEst_E2(fwd_inv_obj, tau, Holo, save_every_f,...
                results_dir_iter,'plotrecon',true,'verbose',true,...
                'x',f_gt,'xhat0',real(BPField));%
            close all;
            % print time in minutes
            toc/60
        end
    end
end

