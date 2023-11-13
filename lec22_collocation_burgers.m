% sample run: lec22_collocation_burgers(100,1.25)

function error=lec22_collocation_wave(N,time_final)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: initialize the grid and initialize u^0_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the grid
x_left = 0.0;
x_right = 2*pi;
h = (x_right-x_left)/N;
x = h*(1:N)';
dx = 1/N;

u0 = sin(x);

dt = dx;
Nt = ceil(time_final/dt);
dt = time_final/Nt;

u_old = u0;
u = u0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: time evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
close all
set(gca,'FontSize',15);
k = [0:N/2-1 0 -N/2+1:-1]';
for time_step = 1:Nt
    flux = 0.5*(u_old.^2);
    flux_hat = fft(flux);
    flux_hat_x = 1i*(k.*flux_hat);
    u = u_old - dt*real(ifft(flux_hat_x));

    u_old = u;

    plot(x,u,'o');
    ylim([-1.5,1.5]);
    xlim([0,2*pi]);
    pause(0.01)
end


end

