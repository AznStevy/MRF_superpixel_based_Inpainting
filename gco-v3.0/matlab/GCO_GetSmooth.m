function u = GCO_GetSmooth(Handle, adj, alpha)
maxiter = 500; epsilon = [100 1]; lambda  = 10; dt = 1;

[M,N,C] = size(adj);
u = adj;

h1 = 1; h2 = 1;
Swap = round(maxiter/2);
ep = [epsilon(1)*ones(numel(1:Swap-1),1); epsilon(2)*ones(numel(Swap:maxiter),1)];
c1 = 1/epsilon(2);
lambda = lambda*alpha;

% Diagonalize the Laplace Operator by: Lu + uL => D QuQ + QuQ D, where 
% Q is nonsingular, the matrix of eigenvectors of L and D is a diagonal matrix.
% We have to compute QuQ. This we can do in a fast way by using the fft-transform:

Lambda1 = spdiags(2*(cos(2*(0:M-1)'*pi/M)-1),0,M,M)/h1^2;
Lambda2 = spdiags(2*(cos(2*(0:N-1)'*pi/N)-1),0,N,N)/h2^2;

Denominator = Lambda1*ones(M,N) + ones(M,N)*Lambda2;

for c = 1:C
    u_hat      = fft2(u(:,:,c));
    lu0_hat    = fft2(lambda(:,:,c).*u(:,:,c));
    
    for it = 1:maxiter
        
        lu_hat     = fft2(lambda(:,:,c).*u(:,:,c));
        Fprime_hat = fft2(2*(2*u(:,:,c).^3-3*u(:,:,c).^2+u(:,:,c)));
        
        u_hat = (dt*(1+lambda(:,:,c)-Denominator/epsilon(2)).*u_hat...
            + dt/ep(it)*Denominator.*Fprime_hat...
            + dt*(lu0_hat-lu_hat))./(1+lambda(:,:,c)*dt+ep(it)*dt*Denominator.^2-dt*Denominator/epsilon(2));
        
        u(:,:,c) = real(ifft2(u_hat));
    end
end
