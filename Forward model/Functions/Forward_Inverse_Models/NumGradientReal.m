function numgrad=NumGradientReal(holo,f,A,eps,numpts)

[Nx,Ny,Nz] = size(f);
numgrad = zeros(Nx,Ny,Nz);

N = numpts;

for iii = 1:N
    % convert linear index to 3D index
    siz = size(f);
    [ii,jj,kk] = ind2sub(siz,iii);
    
    % ei is unit vecor
    ei = zeros(Nx,Ny,Nz);
    ei(ii,jj,kk) = 1;
    
    % this is the wiggle to compute gradient
    f_plus = f + eps * ei;
    f_minus = f - eps * ei;
    
    % J+
    resid = A(f_plus) - holo;
    J_plus = 0.5*norm(resid(:),2)^2;
    
    % J-
    resid = A(f_minus) - holo;
    J_minus = 0.5*norm(resid(:),2)^2;
    
    % copmtue numerical gradient at f(iii)
    grad_re = (J_plus - J_minus) / (2*eps);
    
    %% compute total gradient
    numgrad(iii) = grad_re;
    
    fprintf('iter:%d\n',iii);
end

end