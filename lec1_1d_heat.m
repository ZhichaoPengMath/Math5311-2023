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
% r: dt/dx^2
% T: final time

% run 
% lec1_1d_heat(10,0.5,1.0);
% lec1_1d_heat(20,0.5,1.0);
% lec1_1d_heat(100,0.5,0.5);

% lec1_1d_heat(100,0.51,0.1)
% lec1_1d_heat(100,0.501,0.1);
% lec1_1d_heat(100,0.501,1);

function error_inf=lec_1d_heat(N,r0,T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: initialize the grid and initialize u^0_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the grid
x_left = 0.0;
x_right = 1.0;
x = linspace(0,1,N+1);
dx = 1/N;

% compute dt= r*dx^2
dt = r0*dx^2;
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
% Step 2: time marching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
for time_step = 1:N_time_step
    u_old = u; % get u_j^n

    % Loop over inner points
    for j = 2:N
        u(j) = (1-2*r)*u_old(j)+r*u_old(j-1)+r*u_old(j+1);
    end
    
    % Vectorize version
    %u(2:N) = (1-2*r)*u_old(2:N)+r*u_old(1:N-1)+r*u_old(3:N+1);

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
