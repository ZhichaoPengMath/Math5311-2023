demo_option = 2;

N = 400;
x = linspace(-pi,3*pi,N+1);

dx = 4*pi/N;
dt = dx;
T = pi;
N_time_step = ceil(T/dt);
dt = T/N_time_step;
% Demo of phase speed
close all;
if (demo_option==1)
    for time_step = 1:N_time_step
        time = dt*time_step;
        figure(1)
        plot(x,wave1(x-time),'b-','LineWidth',1.5);
        ylim([-1,1])
        axis tight
        pause(0.05)
    end
elseif (demo_option==2)
    %time = dt*time_step;
    time = 0;
    figure(1)
    h1 = plot(x,wave2(x,time),'b-','LineWidth',1.5);
    hold on
    h2 = plot(x,cos(x-3*time),'r--','LineWidth',1.5);
    xlim([x(1),x(end)]);
    ylim([-1,1]);
    for time_step = 2:N_time_step
        time = dt*time_step;
        figure(1)
        delete(h1)
        delete(h2)
        h1 = plot(x,wave2(x,time),'b-','LineWidth',1.5);
        hold on
        h2 = plot(x,cos(x-3*time),'r--','LineWidth',1.5);
        pause(0.05)
    end
end

function res=wave1(x)
    res = zeros(size(x));
    for i = 1:length(x)
        if (x(i)<2*pi)&&(x(i)>0)
            res(i) = cos(10*x(i));
        end
    end
end

function res = wave2(x,time)
    res = cos(10*(x-5*time)).*cos(x-3*time);
end