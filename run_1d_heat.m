N = [10,20,40,80,160,320];

err_ex = zeros(length(N),1);
err_im = zeros(length(N),1);
err_cn = zeros(length(N),1);

for i = 1:length(N)
    err_ex(i) = lec1_1d_heat_explicit(N(i),0.5,0.5);
    err_im(i) = lec2_1d_heat_implicit(N(i),0.5,0.5);
    err_cn(i) = lec3_1d_heat_crank_nicolson(N(i),0.5,0.5);
end

close all;
figure(1)
loglog(N,err_ex,'-o','Linewidth',1.5);
hold on
loglog(N,0.1*N.^-2,'--','Linewidth',1.5);
legend('Max error','Slope=2');
set(gca,'FontSize',18);
xlabel('N')
xticks(N)
ylim([5e-7,1e-3]);
ylabel('Error')
grid on
title('Explicit')


figure(2)
loglog(N.^-1,err_ex,'-o','Linewidth',1.5);
hold on
loglog(N.^-1,0.1*N.^-2,'--','Linewidth',1.5);
lgd=legend('Max error','Slope=2');
lgd.Location='best';
set(gca,'FontSize',18);
xlabel('dx')
grid on
xticks([1/200,1/100,1/50,1/20,1/10])
ylim([5e-7,1e-3]);
ylabel('Error')
title('Explicit')

figure(3)
loglog(N,err_im,'-o','Linewidth',1.5);
hold on
loglog(N,0.1*N.^-1,'--','Linewidth',1.5);
legend('Max error','Slope=1');
set(gca,'FontSize',18);
xlabel('N')
xticks(N)
%ylim([5e-7,1e-3]);
ylabel('Error')
grid on
title('Implicit')

figure(4)
loglog(N.^-1,err_im,'-o','Linewidth',1.5);
hold on
loglog(N.^-1,0.1*N.^-1,'--','Linewidth',1.5);
lgd=legend('Max error','Slope=1');
lgd.Location='best';
set(gca,'FontSize',18);
xlabel('dx')
grid on
xticks([1/200,1/100,1/50,1/20,1/10])
%ylim([5e-7,1e-3]);
ylabel('Error')
title('Implicit')

figure(5)
loglog(N,err_cn,'-o','Linewidth',1.5);
hold on
loglog(N,0.1*N.^-2,'--','Linewidth',1.5);
legend('Max error','Slope=2');
set(gca,'FontSize',18);
xlabel('N')
xticks(N)
ylim([2e-7,1e-3]);
ylabel('Error')
grid on
title('Crank-Nicolson')

figure(6)
loglog(N.^-1,err_cn,'-o','Linewidth',1.5);
hold on
loglog(N.^-1,0.1*N.^-2,'--','Linewidth',1.5);
lgd=legend('Max error','Slope=2');
lgd.Location='best';
set(gca,'FontSize',18);
xlabel('dx')
grid on
xticks([1/200,1/100,1/50,1/20,1/10])
ylim([2e-7,1e-3]);
ylabel('Error')
title('Crank-Nicolson')

figure(7)
loglog(N,err_cn,'-o','Linewidth',1.5);
hold on
loglog(N,err_ex,'-o','Linewidth',1.5);
loglog(N,0.1*N.^-2,'--','Linewidth',1.5);
legend('CN','ex','Slope=2');
set(gca,'FontSize',18);
xlabel('N')
xticks(N)
ylim([2e-7,1e-3]);
ylabel('Error')
grid on
title('Crank-Nicolson v.s. Explicit')

