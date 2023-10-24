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
% lec3_1d_heat_crank_nicolson(10,0.5,1.0);
% lec3_1d_heat_crank_nicolson(50,0.5,1.0);
% lec3_1d_heat_crank_nicolson(100,0.5,1.0);



function error_inf=lec3_1d_heat_crank_nicolson(N,sigma0,T)
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
% Step 2: Assemble spatial operator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dt/dx^2*()
DiscreteLaplace = zeros(N-1,N-1);
for i = 1:N-1
    DiscreteLaplace(i,i) = -2;
    if (i>1)
        DiscreteLaplace(i,i-1) = 1;
    end

    if(i<N-1)
        DiscreteLaplace(i,i+1) = 1;
    end
end
% Make A a sparse matrix
DiscreteLaplace = sparse(DiscreteLaplace);

Identity = eye(N-1);

% (1-0.5*dt*\partial_{xx} )u_j^{n+1} = (1+0.5*dt*\partial_{xx})u_j^{n}
ImplicitMat = Identity-0.5*r*DiscreteLaplace;
ExplicitMat = Identity+0.5*r*DiscreteLaplace;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: time marching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
for time_step = 1:N_time_step
    u_old = u; % get u_j^n

    rhs_vec = ExplicitMat*u_old(2:N);
    u(2:N) = ImplicitMat\rhs_vec;

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
