% sample run: lec17_jacobi(100);

function error=lec17_jacobi(N)
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
%A = sparse(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Consider a random vector as the exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_exact = rand(N,1);
b = A*u_exact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Jacobi iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_iter = 50000;
residual_history  = zeros(N_iter,1);
residual_history2 = zeros(N_iter,1);
tol = 1e-10;

u_old = zeros(N,1);
u     = zeros(N,1);

% Matrix free implementation
u_old(:) = 0.0;
tic;
for iter = 1:N_iter
    for j = 1:N
        if (j==1)
            u(j) = 0.5*(b(j)+u_old(j+1));
        elseif (j==N)
            u(j) = 0.5*(b(j)+u_old(j-1));
        else
            u(j) = 0.5*(b(j)+u_old(j+1)+u_old(j-1));
        end
    end
    
    u_old = u;
    res = A*u-b;
    if (norm(res)/norm(b)<1e-10)
        break;
    end
    residual_history(iter) = norm(res)/norm(b);
end
toc;
error = norm(u_exact-u,'inf');

% An implementation using Dense matrix
u_old(:) = 0.0;
D = diag(diag(A));
LplusU = A-D;

tic;
for iter = 1:N_iter
    u = D\(b-LplusU *u_old);
    u_old = u;
    res = A*u-b;
    if (norm(res)/norm(b)<1e-10)
        break;
    end
    residual_history(iter) = norm(res)/norm(b);
end
toc;
error2 = norm(u_exact-u,'inf');

close all
figure(1)
semilogy(residual_history(1:iter),'-o');

fprintf('iter = %d, error = %e\n,',iter,error);

end