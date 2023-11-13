function error=lec21_cg_smoothing(N,N_iter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: initialize the grid and initialize u^0_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the grid
x_left = 0.0;
x_right = 1.0;
x = linspace(x_left,x_right,N+1);
% x is a row vector, we want u to be a column vetor so take transpose here
x = x';
x_interior = x(2:end-1);

dx = 1/N;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Assemble spatial discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(N,N);
for j = 1:N
    if (j>1)
        A(j,j-1) = -1;
    end

    if (j<N)
        A(j,j+1) = -1;
    end

    A(j,j) = 2;
end
A = sparse(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Consider a random vector as the exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(10086)
u_exact = rand(N,1);
b = A*u_exact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Jacobi iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
residual_history  = zeros(N_iter,1);
residual_history2 = zeros(N_iter,1);
tol = 1e-10;

u = zeros(N,1);

% Matrix free implementation
u_old(:) = 0.0;
close all
figure(1)

tic;
tmp = A*u;
r = b-tmp;
p = r;
rtr = dot(r,r);
rtr0 = rtr;
for iter = 1:N_iter
    plot(u-u_exact,'-bo','Linewidth',1.5)
    q = A*p;
    alpha = rtr/dot(p,q);
    u = u + alpha*p;
    r = r - alpha*q;
    rtrold = rtr;
    rtr = dot(r,r);

    beta = rtr/rtrold;
    p = beta*p + r;


    if (mod(iter,5)==0)
       plot(u-u_exact,'-bo','Linewidth',1.5)
    end
    pause(0.01)

    if (norm(r)/norm(b) < tol)
        break;
    end

end
toc;
error = norm(u_exact-u,'inf');
fprintf('iter = %d, error = %e\n,',iter,error);

end