% Solve u_t + u_x = 0 with a square pulse

% Run lec8_nonconservative(600,0.5,1);


% flux_option = 1: Godunov
% flux_option = 2: Central
function lec8_nonconservative(N,sigma,T)
global dx
global dt
% Step 1: Define the mesh
x_left = 0;
x_right = 3;
x_edge = linspace(x_left,x_right,N+1);
x = 0.5*( x_edge(1:end-1)+x_edge(2:end) );
% I prefer a column vector
x = x';
dx = (x_right-x_left)/N;
dt = sigma*dx;
N_time_step = ceil(T/dt);
dt = T/N_time_step;
sigma = dt/dx;

% Step 2: initialize the scheme
u = u_initial(x);
u_old = u;

% Step 3: loop over time
close all
figure(1)
h1 = plot(x,u_initial(x),'r-','Linewidth',1.5);
hold on 
h2 = plot(x,u_old,'bo');
    ylim([-0.5,1.5]);
    xlim([x_left,x_right]);


for time_step = 1:N_time_step
    u_old = u;
    % interior update, zero boundary conditions are exactly imposed
    % flux for the left boundary
    flux_right = 0; 
    for i = 1:N
        if (i==1)
            u_left = 0.0;
        else
            u_left = u_old(i-1);
        end
        u(i) = u_old(i) - dt/dx*u_old(i)*(u_old(i)-u_left);

    end

    % Plot
    delete(h1)
    delete(h2)
    h1 = plot(x,u_exact(x,time_step*dt),'r-','Linewidth',1.5);
    hold on 
    h2 = plot(x,u,'bo');
    pause(0.01);
end
%%%%%
end

function res = u_initial(x)
    res = zeros(size(x));
    for i = 1:length(x)
        if (x(i)>1)&&(x(i)<2)
            res(i) = 1.0;
        end
    end
end

function res = u_exact(x,time)
    res = zeros(size(x));
    for i = 1:length(x)
        if ( x(i)<=2+0.5*time)&&(x(i)>=1+time)
           res(i) = 1.0; 
        elseif ( x(i)>=1 )&&(x(i)<1+time)
            res(i) = (x(i)-1)/time;
        end
    end
end



function res = numerical_flux(u_l,u_r,flux_option)
    if (flux_option==1)
        if (u_l<u_r)
            if (u_l>0)||(u_r<0)
                res = min(0.5*u_l^2,0.5*u_r^2);
            else
                res = 0.0;
            end
        else
            res = max(0.5*u_l^2,0.5*u_r^2);
        end
    else
        res = 0.25*(u_l+u_r)^2;
    end
end
