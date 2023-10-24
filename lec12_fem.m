%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo code for solving:
%      -u_xx=2, x in (0,1)
%  IC: u(x,0) = sin(pi*x),
%  BC: u(0,t) = u(1,t)=0.
%
% Exact solution to this problem:
%    u(x,t) = sin(pi*x)*exp(-pi^2*t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N: Number of interior points
% N+1: Number of elements


function error_inf=lec12_fem(N,option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: initialize the grid and initialize u^0_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the grid
x_left = 0.0;
x_right = 1.0;
x = linspace(x_left,x_right,N+2);
% x is a row vector, we want u to be a column vetor so take transpose here
x = x';
x_interior = x(2:end-1);

dx = zeros(N+1,1);
for i = 1:N+1
    dx(i) = x(i+1)-x(i);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Assemble spatial discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(N,N);
for j = 1:N
    if (j>1)
        A(j,j-1) = -1/dx(j);
    end

    if (j<N)
        A(j,j+1) = -1/dx(j+1);
    end

    A(j,j) = 1/dx(j)+1/dx(j+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Assemble right hand side, -u_xx = 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = zeros(N,1);
for j = 1:N
    f(j) = source_function(x_interior(j),option)*0.5*(dx(j)+dx(j+1));

    % Mid point
    %f(j) = 0.5*(source_function(x_interior(j)-0.5*dx(j),option)*dx(j)...
    %           +source_function(x_interior(j)+0.5*dx(j),option)*dx(j+1));

    % Linear interpolation, note that  x_interior(j) = x(j+1) 
    %{
    f(j) = 1/6*(   dx(j)*source_function(x(j),option)...
               +2*(dx(j)+dx(j+1))*source_function(x(j+1),option)...
                  +dx(j+1)*source_function(x(j+2),option) );
    %}
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = A\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: demonstrate the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute exact solution
u_ex = u_exact(x_interior,option);
% compare the exact solution and the numerical solution
error_inf = max(abs(u_ex-u));
fprintf("Max error is %e \n",error_inf);

close all;
plot(x_interior,u_ex,'r-','Linewidth',1.5);
hold on
plot(x_interior,u,'bo','Linewidth',1.5);
legend('Exact','Numerical solution')
xlabel('x');
ylabel('y');
font_size = 15;
set(gca,'FontSize',font_size);
box on

end

% source function
function res = source_function(x,option)
    if (option == 1)
        res = 2;
    elseif (option==2)
        res = pi^2*sin(pi*x);
    end
end

% exact solution
function res = u_exact(x,option)
    if (option == 1)
        res = x.*(1-x);
    else
        res = sin(pi*x); 
    end
end
