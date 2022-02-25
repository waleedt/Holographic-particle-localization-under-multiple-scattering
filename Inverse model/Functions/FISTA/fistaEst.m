function [xhat, outs] = fistaEst(forwardObj, tau, Holo2D, fileID, results_dir, varargin)
% FISTAEST for solving TV optimization problems.
%
% Input:
%   - forwardObj: forward model operator
%   - tau/sv_sq: regularization parameter
%   - Holo2D: 2D measured hologram (after DC subraction)
%   - warmup_iters: #iterations to run FPPA even if cost
%     is initially increasing
%   - varargin: options for the algorithm
%       Options:
%       - numIter: number of iterations (def: 100)
%       - plotRecon: plot reconstruction (def: false)
%       - step: step-size of the algorithm
%               (recommendation is to use 1/Lipschitz constant of H)
%       - tol: stopping criterion for the algorithm (def: 1e-4)
%       - tolCount: number of satisfactions of tol before stop (def: 10)
%       - verbose: print command line message (def: false)
%       - x: oracle i.e. the original signal
%       - xhat0: initial value
%
% Output:
% - xhat: estimated signal
% - outs: extra data
%       - outs.cost: per iteration cost
%       - outs.snr: per iteration snr
%       - outs.time: per iteration CPU-time
%       - outs.relativeChange: per iteration distance from tolerance
%
% U. S. Kamilov, MERL, 2015.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigSize = forwardObj.sigSize; % extract the signal size

eta = 1/2;
backTol = 1e-12;
numIter = 100;
plotRecon = true;
step = forwardObj.step;
tol = 1e-4;             %waleed: set appropriate tol for stop criterion
cost_tol = 1e-6;
tolCount = 10;
cost_tolCount = 10;
verbose = true;
xhat0 = zeros(sigSize); %waleed: do initialization with At(x)
save_f_every_x_iters = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nargs = length(varargin); % Number of options

for i = 1:2:nargs % Go through options
    
    name = lower(varargin{i}); % Extract name
    value = varargin{i+1}; % Extract value
    
    switch(name)
        case 'numiter'
            numIter = value;
        case 'plotrecon'
            plotRecon = value;
        case 'step'
            step = value;
        case 'tol'
            tol = value;
        case 'costtol'
            cost_tol = value;
        case 'tolcount'
            tolCount = value;
        case 'verbose'
            verbose = value;
        case 'x'
            x = value;
        case 'xhat0'
            xhat0 = value;
        otherwise
            error('fistaEst: input is not recognized!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isOracleSet = exist('x', 'var'); % true signal is set

%%% Create gradient object for TV
gradObj = GradientOperator3D(sigSize);

%%% useful handles
tvCost = @(x) sum(sum(sum(sqrt(sum(abs(gradObj.mult(x)).^2, 4)))));

evaluateSnr = @(x, xhat) 20*log10(norm(x(:))/norm(x(:)-xhat(:)));
evaluateTol = @(x,xnext) norm(x(:)-xnext(:))/norm(x(:));

%%% Initialize
xhat = forwardObj.AT_init(Holo2D);

% enforce positivity on initial (least sq) solution
% if real/imag(xhat)<0 --> real/imag(xhat)=0
%%%xhat(imag(xhat)< 0) = xhat(imag(xhat)< 0) - 1i*imag(xhat(imag(xhat)< 0));
%%%xhat(real(xhat)< 0) = xhat(real(xhat)< 0) - real(xhat(real(xhat)< 0));

xhat(:) = 0;

s = xhat;
q = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track diagnostic information [time, cost, SNR]
if(plotRecon)
    outs.time = zeros(numIter, 1);
    outs.cost = zeros(numIter, 1);
    
    h = figure('Name', 'fistaEst');
end

if(isOracleSet)
    outs.snr = zeros(numIter, 1);
    snr_all = zeros(50000,1);
    snr_all(1) = evaluateSnr(x, xhat);
    
    outs.relativeChange = zeros(numIter, 1);
    relativeChange_all = zeros(50000,1);
    relativeChange_all(1) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tolCounter = tolCount; % conter for tolerance violations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[z_est,SF,E_rb_all,~] = forwardObj.A(s);
resid_2D = z_est - Holo2D;
[cost_init,fid_cost_init,reg_cost_init] = evaluateCost2(resid_2D,xhat,tau,tvCost);
iIter1 = 1;
outs.cost(iIter1) = cost_init;

cost_all = zeros(50000,1);
cost_fid = zeros(50000,1);
cost_reg = zeros(50000,1);
cost_all(iIter1) = cost_init;
cost_fid(iIter1) = fid_cost_init;
cost_reg(iIter1) = reg_cost_init;

gradNorm = zeros(50000,1);
relativeChange_grad = zeros(50000,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plotRecon)
    figure(h);
    set(h, 'Name', sprintf('fistaEst [%d/%d]', iIter1, numIter));
    
    semilogy(1:iIter1, outs.cost(1:iIter1), 'b-',...
        iIter1, outs.cost(iIter1), 'ro', 'LineWidth', 1.5);
    xlim([1 numIter]);
    grid on;
    hold on;
    title(sprintf('cost: %.4e', outs.cost(iIter1)));
    set(gca, 'FontSize', 14);
    drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cost_prev = cost_init;
grad_prev = zeros(size(xhat));
msg = sprintf('Initial Cost Value = %2.10e\n',cost_init);
% if(verbose)
%     fprintf(msg);
% end
fprintf(fileID,msg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iIter = 1;
while(1)
    tic
    timeStart = cputime;
    
    %%% proximal-Gradient Step
    grad = forwardObj.Grad(resid_2D,s,E_rb_all);
    xhatnext = proxTV(s-step*grad, step*tau, gradObj);
    xhatnext(xhatnext<0) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All additional code for backtracking goes here (waleed)
    d = evaluateFidCost(s,Holo2D,forwardObj);
    [data,z_est,SF,E_rb_all] = evaluateFidCost_BTR(xhatnext,Holo2D,forwardObj);
    quad = d...
        + (0.5/step)*(norm(xhatnext(:)-s(:)+step*grad(:))^2)...
        - (0.5*step)*(norm(grad(:))^2);
    
    % Backtracking check the line search condition
    while(backTol > 0 && (quad - data) < 0)
        % reduce step
        step = eta*step;
        
        % tentative update
        xhatnext = proxTV(s-step*grad, step*tau, gradObj);
        
        % data-term and quadratic upper bound
        [data,z_est,SF,E_rb_all] = evaluateFidCost_BTR(xhatnext,Holo2D,forwardObj);
        quad = d...
            + (0.5/step)*(norm(xhatnext(:)-s(:)+step*grad(:))^2)...
            - (0.5*step)*(norm(grad(:))^2);
        
        % print to command line
        if(verbose)
            fprintf('[fistaIter: %d]', iIter);
            fprintf('[update attempt]');
            fprintf('[step: %.1e]', step);
            fprintf('[data: %.1e]', data);
            fprintf('[quad: %.1e]', quad);
            fprintf('[diff: %.1e]', quad-data);
            fprintf('\n');
        end
        
        % if step is smaller than tolerance
        if(step < backTol)
            warning('fistaEst: iter %d of %d: step < %.1e',...
                indIter, numIter, backTol);
            step = backTol;
            break;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% FISTA over-relaxation
    qnext = 0.5*(1+sqrt(1+4*q*q));
    s = xhatnext + ((q-1)/qnext)*(xhatnext-xhat);
    
    % Update residual for cost and grad computation
    %[z_est,SF,E_rb_all,~] = forwardObj.A(xhatnext);
    resid_2D = z_est - Holo2D;
    
    %%% Update diagnistic information
    fid_cost_now = data;
    [cost_now,reg_cost_now] = evaluateCost3(fid_cost_now,xhatnext,tau,tvCost);
    
    if(plotRecon)
        outs.cost(iIter+1) = cost_now;
    end
    
    cost_all(iIter+1) = cost_now;
    cost_fid(iIter+1) = fid_cost_now;
    cost_reg(iIter+1) = reg_cost_now;
    
    relativeChange_grad(iIter) = evaluateTol(grad, grad_prev);
    gradNorm(iIter) = norm(grad(:));
    grad_prev = grad;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(plotRecon)
        outs.time(iIter) = cputime - timeStart;
    end
    
    if(isOracleSet)
        outs.snr(iIter) = evaluateSnr(x, xhatnext);
        snr_all(iIter+1) = evaluateSnr(x, xhatnext);
        outs.relativeChange(iIter) = evaluateTol(xhat, xhatnext);
        relativeChange_all(iIter+1) = evaluateTol(xhat, xhatnext);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % update log file and disp progress
    msg = sprintf('***FISTA Iteration %d***\n',iIter);
    fprintf(fileID,msg);
    
    %     if(verbose)
    %         fprintf(msg);
    %         fprintf('Iter: %d , Iter time: %f seconds\n',iIter,toc);
    %         fprintf('Born: %d, Tau: %d\n',forwardObj.bornorder,tau);
    %     end
    
    if (verbose)
        fprintf('[fistaIter: %d]', iIter);
        if(cost_now < cost_prev)
            msg = sprintf('[reduced cost: %2.10e]\n',cost_now);
        elseif (cost_now > cost_prev)
            msg = sprintf('[increased cost: %2.10e]\n',cost_now);
        else
            msg = sprintf('[equal cost: %2.10e]\n',cost_now);
        end
        fprintf(msg);
    end
    
    fprintf(fileID,msg);
    
    
    % update objective cost
    delta_cost = cost_now - cost_prev;
    cost_prev = cost_now;
    if(~isfinite(cost_now))
        msg = sprintf('cost -> infinite. Exit\n');
        fprintf(fileID,msg);
        break;
    end
    
    %%% Update previous iterates
    q = qnext;
    xhat = xhatnext;
    
    %%% Print diagnostics (Not used here)
    %    if(verbose)
    %         fprintf('[fistaEst: %d/%d]', iIter, numIter);
    %         fprintf('[cost: %.2e]', outs.cost(iIter+1));
    %         fprintf('[time: %.1f]', sum(outs.time));
    %         fprintf('[tols: %.4f]', outs.relativeChange(iIter));
    %         if(isOracleSet)
    %             fprintf('[snr: %.2f]', outs.snr(iIter));
    %         end
    %         fprintf('\n');
    %        fprintf('[cost: %.2e]\n', cost_all(iIter1));
    
    %    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot diagnostics
    if(plotRecon)
        figure(h);
        set(h, 'Name', sprintf('fistaEst [%d/%d]', iIter+1, numIter));
        
        if(~isOracleSet)
            semilogy(1:iIter, outs.cost(1:iIter), 'b-',...
                iIter+1, outs.cost(iIter+1), 'ro', 'LineWidth', 1.5);
            xlim([1 numIter]);
            grid on;
            title(sprintf('cost: %.4e', outs.cost(iIter+1)));
            set(gca, 'FontSize', 14);
        else
            subplot(1, 2, 1);
            semilogy(1:iIter, outs.cost(1:iIter), 'b-',...
                iIter, outs.cost(iIter), 'ro', 'LineWidth', 1.5);
            xlim([1 numIter]);
            grid on;
            title(sprintf('cost: %.4e', outs.cost(iIter)));
            set(gca, 'FontSize', 14);
            
            subplot(1, 2, 2);
            plot(1:iIter, outs.snr(1:iIter), 'b-',...
                iIter, outs.snr(iIter), 'ro', 'LineWidth', 1.5);
            xlim([1 numIter]);
            grid on;
            title(sprintf('snr: %.2f', outs.snr(iIter)));
            set(gca, 'FontSize', 14);
        end
        drawnow;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save diagnostics and break
    %     if(tolCounter <= 0)
    %         if(isOracleSet)
    %             outs.snr = outs.snr(1:iIter);
    %         end
    %         outs.cost = outs.cost(1:iIter);
    %         outs.relativeChange= outs.relativeChange(1:iIter);
    %         outs.time = outs.time(1:iIter);
    %         msg = sprintf('tolCounter <= 0.Exit.\n');
    %         fprintf(fileID,msg);
    %         break;
    %     end
    
    % save object data after every 100 iterations
    if(rem(iIter,save_f_every_x_iters)==0)
        save_intermediate_f;
    end
    
%     %% Stopping criterion relative change
%         if(outs.relativeChange(iIter) < tol)
%             tolCounter = tolCounter - 1;
%         else
%             tolCounter = tolCount;
%         end
%         if(tolCounter == 0)
%             if(verbose)
%                 fprintf('delta(cost) < cost_tol.Exit.\n');
%             end
%             break;
%         end

    %%% Stopping criterion delta_cost
    if(abs(delta_cost) <= cost_tol)
        cost_tolCount = cost_tolCount - 1;
        if(cost_tolCount == 0)
            msg = sprintf('delta(cost) < cost_tol.Exit.\n');
            fprintf(fileID,msg);
            break;
        end
    end
    toc
    iIter = iIter + 1;
end

% save cost vs iteration figure
% if(plotRecon)
%     born_order = forwardObj.bornorder;
%     mkdir(sprintf('Results/born%d_cost',born_order));
%     mkdir(sprintf('%s/cost',results_dir));
%     saveas(h,sprintf('Results/born%d_cost/cost_tau_%1.1e_b%d.png',born_order,tau,born_order),'png');
%     saveas(h,sprintf('%s/cost/cost_tau_%1.1e_b%d.png',results_dir,tau,born_order),'png');
% end
save_intermediate_f;
end