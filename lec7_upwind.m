% Solve u_t + u_x = 0 with a square pulse
% Run lec7_upwind(100,0.75,2.5);
% Run lec7_upwind(1000,0.75,2.5);

function lec7_upwind(N,sigma,T)
% Step 1: Define the mesh
x_left = 0;
x_right = 10;
x = linspace(x_left,x_right,N+1);
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
h2 = plot(x,u,'bo');
ylim([-0.5,1.5]);
xlim([0,10]);
for time_step = 1:N_time_step
    u_old = u;
    % interior update, zero boundary conditions are exactly imposed
    ind_interior = 2:N;
    u(ind_interior) = u_old(ind_interior)...
              -sigma*(u_old(ind_interior)-u_old(ind_interior-1));

    % Plot
    delete(h1)
    delete(h2)
    h1 = plot(x,u_initial(x-time_step*dt),'r-','Linewidth',1.5);
    hold on 
    h2 = plot(x,u,'bo');
    pause(0.05);
end


end

function res = u_initial(x)
    res = zeros(size(x));
    for i = 1:length(res)
        if (x(i)>1)&&(x(i)<2)
            res(i) = 1.0;
        end
    end
end
