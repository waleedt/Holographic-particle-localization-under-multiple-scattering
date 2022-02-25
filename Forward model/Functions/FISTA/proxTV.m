%xhatnext             = proxTV(xhatnext, step*tau , gradObj);
function [xhat, outs] = proxTV(y       , lam      , gradObj)

%%% TV denoising of a complex image 
%%%
%%% U. S. Kamilov

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

computeSnr = @(x, xhat) 20*log10(norm(MyV2C(x))/norm(MyV2C(x)-MyV2C(xhat)));
tvCost = @(x) sum(sum(sum(sqrt(sum(abs(gradObj.mult(x)).^2, 4)))));
computeCost = @(x) 0.5*(norm(MyV2C(x)-MyV2C(y))^2)+lam*tvCost(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numIter = 10; % Number of iterations
plotRecon = false; % Plot the progress
tol = 1e-6; % Tolerance for relative change
verbose = false; % Display messages

rho = 1; % quadratic parameter
xhat0 = y; % initial image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial values
xhat = xhat0;
d = gradObj.mult(xhat);
xhatmult = d;
s = zeros(size(d));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterative reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iIter = 1:numIter
    
    % message to command line
    if(verbose)
        fprintf('[proxTV: %d/%d]', iIter, numIter);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update estimates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xhatprev = xhat; % store the current estimate before update
    
    % update d
    d = shrink3D_real(xhatmult-s/rho, lam/rho, 0);

    % update xhat
    dat = y + rho*gradObj.multTranspose(d + s/rho);
    xhat = real(ifftn(fftn(dat) ./ (1 + rho*gradObj.freqResponseDTD)));
    xhatmult = gradObj.mult(xhat);
    % update s
    s = s + rho*(d - xhatmult);
    
    % check tolerance
    if(norm(xhat(:)-xhatprev(:))/norm(xhatprev(:)) < tol)
        break;
    end
end