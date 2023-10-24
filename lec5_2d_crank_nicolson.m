%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo code for solving:
%      u_t = u_xx+u_yy, (x,y) in (0,pi)^2
%  IC: u(x,0) = sin(x)*sin(y),
%  BC: u(0,t) = u(pi,t)=0.
%
% Exact solution to this problem:
%    u(x,t) = sin(x)*sin(y)*exp(-2*t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N: how many grid point, we use in space
% sigma0: dt/dx
% T: final time

function error_inf=lec5_2d_crank_nicolson(Nx,Ny,sigma0,T)

assemble_option=0;
lu_option = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: initialize the grid and initialize u^0_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the grid
x_left = 0.0;
x_right = pi;
y_left = 0.0;
y_right = pi;
x = linspace(x_left,x_right,Nx+1);
y = linspace(y_left,y_right,Ny+1);
dx = (x_right-x_left)/Nx;
dy = (y_right-y_left)/Ny;
h = min(dx,dy);

% compute dt= r*dx^2
dt = sigma0*h;
% Make sure N_time_step*dt = T
N_time_step = ceil(T/dt);
dt = T/N_time_step;

% Get rx,ry
rx = dt/dx^2;
ry = dt/dy^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: define a matrix to save index map and initialize
% the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_map = zeros(Nx-1,Ny-1);
ind = 1;
for j = 1:Nx-1
    for i = 1:Ny-1
        ind_map(i,j) = ind;
        ind = ind+1;
    end
end

tot_dof = (Nx-1)*(Ny-1);
u = zeros( tot_dof,1 );
for i = 1:Nx-1
    for j = 1:Ny-1
        u(ind_map(i,j)) = u_initial( x(i+1),y(j+1) );
    end
end
u_old = u;
rhs_vec = u;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: assemble spatial operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (assemble_option==0)
    A_implicit = zeros(tot_dof,tot_dof);
    A_explicit = zeros(tot_dof,tot_dof);

    tic;
    for i = 1:Nx-1
        for j = 1:Ny-1
            A_implicit( ind_map(i,j),ind_map(i,j) ) = 1+rx+ry;
            A_explicit( ind_map(i,j),ind_map(i,j) ) = 1-rx-ry;
        
            if (i>1)
                A_implicit( ind_map(i,j),ind_map(i-1,j) ) = -0.5*rx;
                A_explicit( ind_map(i,j),ind_map(i-1,j) ) =  0.5*rx;
            end

            if (i<Nx-1)
                A_implicit( ind_map(i,j),ind_map(i+1,j) ) = -0.5*rx;
                A_explicit( ind_map(i,j),ind_map(i+1,j) ) =  0.5*rx;
            end

            if (j>1)
                A_implicit( ind_map(i,j),ind_map(i,j-1) ) = -0.5*ry;
                A_explicit( ind_map(i,j),ind_map(i,j-1) ) =  0.5*rx;
            end

            if (j<Ny-1)
                A_implicit( ind_map(i,j),ind_map(i,j+1) ) = -0.5*ry;
                A_explicit( ind_map(i,j),ind_map(i,j+1) ) =  0.5*rx;
            end
        end
    end
    toc;
    A_implicit = sparse(A_implicit);
    A_explicit = sparse(A_explicit);
else
    tic;
    row_ind = zeros(tot_dof*5);
    col_ind = zeros(tot_dof*5);
    val_implicit = zeros(tot_dof*5);
    val_explicit = zeros(tot_dof*5);

    ind = 0;
    for i = 1:Nx-1
        for j = 1:Nx-1
            ind = ind+1;
            row_ind(ind) = ind_map(i,j);
            col_ind(ind) = ind_map(i,j);
            val_implicit(ind) = 1+rx+ry;
            val_explicit(ind) = 1-rx-ry;
        
            if (i>1)
                ind = ind+1;
                row_ind(ind) = ind_map(i,j);
                col_ind(ind) = ind_map(i-1,j);
                val_implicit(ind) = -0.5*rx;
                val_explicit(ind) =  0.5*rx;
            end

            if (i<Nx-1)
                ind = ind+1;
                row_ind(ind) = ind_map(i,j);
                col_ind(ind) = ind_map(i+1,j);
                val_implicit(ind) = -0.5*rx;
                val_explicit(ind) =  0.5*rx;
            end

            if (j>1)
                ind = ind+1;
                row_ind(ind) = ind_map(i,j);
                col_ind(ind) = ind_map(i,j-1);
                val_implicit(ind) = -0.5*ry;
                val_explicit(ind) =  0.5*ry;
            end

            if (j<Ny-1)
                ind = ind+1;
                row_ind(ind) = ind_map(i,j);
                col_ind(ind) = ind_map(i,j+1);
                val_implicit(ind) = -0.5*ry;
                val_explicit(ind) =  0.5*ry;
            end
        end
    end
    A_implicit = sparse( row_ind(1:ind),col_ind(1:ind),val_implicit(1:ind) );
    A_explicit = sparse( row_ind(1:ind),col_ind(1:ind),val_explicit(1:ind) );
    toc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: time marching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
for time_step = 1:N_time_step
    u_old = u; % get u_j^n

    % Loop over inner points to assemble right hand side
    rhs_vec = A_explicit*u_old;
    u = A_implicit\rhs_vec;
    
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
u_ex = u;
for i = 1:Nx-1
    for j = 1:Ny-1
        u_ex(ind_map(i,j)) = u_exact( x(i+1),y(j+1),T );
    end
end
error_inf = max(abs(u_ex-u));
fprintf("Max error is %e \n",error_inf);

[X,Y] = meshgrid(x(2:Nx),y(2:Ny));
u_plot = reshape(u,[Nx-1,Ny-1]);
u_ex_plot = reshape(u_ex,[Nx-1,Ny-1]);
close all;

figure(1)
surf(X,Y,u_plot');
shading interp
view(2)
box on
xlabel('x')
ylabel('y')
set(gca,'FontSize',15);
axis equal;
axis tight;
colorbar;
title('Numerical solution')

figure(2)
surf(X,Y,u_ex_plot');
shading interp
view(2)
box on
xlabel('x')
ylabel('y')
set(gca,'FontSize',15);
axis equal;
axis tight;
colorbar;
title('Exact solution')

figure(3)
surf(X,Y,u_ex_plot');
shading interp
view(2)
box on
xlabel('x')
ylabel('y')
set(gca,'FontSize',15);
axis equal;
axis tight;
colorbar;
title('Error')

end



% initial condition
function res = u_initial(x,y)
    res = sin(x)*sin(y);
end

% exact solution
function res = u_exact(x,y,t)
    res = exp(-2*t)*sin(x)*sin(y);
end
