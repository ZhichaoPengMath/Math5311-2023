%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo code for solving:
%      u_t = u_xx, x in (0,1)
%  IC: u(x,0) = sin(pi*x),
%  BC: u(0,t) = u(1,t)=0.
%
% Exact solution to this problem:
%    u(x,t) = sin(pi*x)*exp(-pi^2*t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N: how many grid point, we use in space
% sigma: dt/dx
% T: final time

% run 
% lec2_1d_heat_implicit(10,0.5,1.0);
% lec2_1d_heat_implicit(20,0.5,1.0);
% lec2_1d_heat_implicit(100,0.5,0.5);

% lec2_1d_heat_implicit(100,0.51,0.1)
% lec2_1d_heat_implicit(100,0.501,0.1);
% lec2_1d_heat_implicit(100,0.501,1);

function error_inf=lec2_1d_heat_implicit(N,sigma0,T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: initialize the grid and initialize u^0_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the grid
x_left = 0.0;
x_right = 1.0;
x = linspace(x_left,x_right,N+1);
% x is a row vector, we want it to be a column vector for convenience
x = x';
dx = (x_right-x_left)/N;

% compute dt= r*dx^2
dt = sigma0*dx;
% Make sure N_time_step*dt = T
N_time_step = ceil(T/dt);
dt = T/N_time_step;
% Get r
r = dt/dx^2;

% allocate memory for u^{n+1} and u^n
% and impose initial conditions
u = u_initial(x); 
u_old = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Assemble linear system 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = zeros(N-1,N-1);
for i = 1:N-1
    A(i,i) = 1+2*r;
    if (i-1>0)
        A(i,i-1) = -r;
    end

    if(i+1<N)
        A(i,i+1) = -r;
    end
end
% Make A a sparse matrix
A = sparse(A);


% Directly use sparse format
%{
row_ind = zeros(3*(N-1),1);
col_ind = zeros(3*(N-1),1);
val_ind = zeros(3*(N-1),1);

ind = 0;
for i = 1:N-1
    ind = ind+1;
    row_ind(ind) = i;
    col_ind(ind) = i;
    val_ind(ind) = 1+2*r;

    if (i-1>0)
        ind = ind+1;
        row_ind(ind) = i;
        col_ind(ind) = i-1;
        val_ind(ind) = -r;
    end

    if (i+1<N-1)
        ind = ind+1;
        row_ind(ind) = i;
        col_ind(ind) = i+1;
        val_ind(ind) = -r;
    end
end
A = sparse(row_ind(1:ind),col_ind(1:ind),val_ind(1:ind));
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: time marching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
for time_step = 1:N_time_step
    u_old = u; % get u_j^n

    rhs_vec = u_old(2:N);
    u(2:N) = A\rhs_vec;

    % impose BC
    % u(N+1)=0;
    % u(N)=0;

end
toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: demonstrate the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute exact solution
u_ex = u_exact(x,T);
% compare the exact solution and the numerical solution
error_inf = max(abs(u_ex-u));
fprintf("Max error is %e \n",error_inf);

close all;
plot(x,u_ex,'r-','Linewidth',1.5);
hold on
plot(x,u,'b--','Linewidth',1.5);
legend('Exact','Numerical solution')
xlabel('x');
ylabel('y');
font_size = 15;
set(gca,'FontSize',font_size);
box on

end

% initial condition
function res = u_initial(x)
    res = sin(pi*x);
end

% exact solution
function res = u_exact(x,t)
    res = exp(-pi^2*t)*sin(pi*x);
end
