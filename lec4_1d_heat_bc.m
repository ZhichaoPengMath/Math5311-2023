%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo code for solving:
%      u_t = u_xx, x in (0,1)
%  IC: u(x,0) = sin(pi*x),
%  BC: u_x(0,t) = u(0,t)+pi*exp(-pi^2*t);
%      u(1,t)=0.
%
% (u_{1}-u_{-1})/2/dx = u_0+g(t)
%  u_{-1} = u_{1}-2*dx*u_{0}-2*dx*g(t)
% Exact solution to this problem:
%    u(x,t) = sin(pi*x)*exp(-pi^2*t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N: how many grid point, we use in space
% sigma: dt/dx
% T: final time

% run 
% lec4_1d_heat_bc(50,0.5,0.5);
% lec4_1d_heat_bc(100,0.5,0.5);



function error_inf=lec4_1d_heat_bc(N,sigma0,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: initialize the grid and initialize u^0_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the grid
x_left =  0.0;
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
% u_{-1} = u_{1}-2*dx*u_0-2*dx*left_bc(t)
ind_interior = 1:N;
DiscreteLaplace = zeros(N,N);
for i = 1:N
    DiscreteLaplace(i,i) = -2;
    if (i>1)
        DiscreteLaplace(i,i-1) = 1;
    end

    if(i<N)
        DiscreteLaplace(i,i+1) = 1;
    end

    % Impose boundary condition
    if (i==1)
        DiscreteLaplace(i,i+1) = DiscreteLaplace(i,i+1)+1;
        DiscreteLaplace(i,i) = DiscreteLaplace(i,i)-2*dx;
    end
end
% Make A a sparse matrix
DiscreteLaplace = sparse(DiscreteLaplace);

Identity = sparse(eye(N));

% (1-0.5*r*\partial_{xx} )u_j^{n+1} = (1+0.5*r*\partial_{xx})u_j^{n}
ImplicitMat = Identity-0.5*r*DiscreteLaplace;
ExplicitMat = Identity+0.5*r*DiscreteLaplace;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: time marching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
time = 0;
close all;
for time_step = 1:N_time_step
    u_old = u; % get u_j^n

    rhs_vec = ExplicitMat*u_old(ind_interior);

    % Impose left bc from 0.5*r*(u_{-1}^n+u_{-1}^{n+1})
    rhs_vec(1) =rhs_vec(1)...
          -r*dx*(g(time)+g(time+dt));%(left_bc(time)+left_bc(time+dt));

    % update u
    u(ind_interior) = ImplicitMat\rhs_vec;
    
    % update time
    time = time+dt;
end
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: demonstrate the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute exact solution
u_ex = u_exact(x,time);
% compare the exact solution and the numerical solution
error_inf = max(abs(u_ex(ind_interior)-u(ind_interior)));
%error_inf = norm(u_ex(ind_interior)-u(ind_interior))*sqrt(dx);
fprintf("Max error is %e \n",error_inf);

close all;
plot(x(ind_interior),u_ex(ind_interior),'r-','Linewidth',1.5);
hold on
plot(x(ind_interior),u(ind_interior),'b--','Linewidth',1.5);
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

% boundary condition
function res = g(t)
    res = pi*exp(-pi^2*t);%-u_exact(-0.5,t);%pi*exp(-pi^2*t);
end


% exact solution
function res = u_exact(x,t)
    res = exp(-pi^2*t)*sin(pi*x);
end
