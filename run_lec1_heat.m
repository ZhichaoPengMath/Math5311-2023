N = [10,20,40,80,160,320];

err = zeros(length(N),1);

for i = 1:length(N)
    err(i) = lec1_1d_heat(N(i),0.5,0.5);
end

close all;
figure(1)
loglog(N,err,'-o','Linewidth',1.5);
hold on
loglog(N,0.1*N.^-2,'--','Linewidth',1.5);
legend('Max error','Slope=2');
set(gca,'FontSize',18);
xlabel('N')
xticks(N)
ylim([5e-7,1e-3]);
ylabel('Error')
grid on

figure(2)
loglog(N.^-1,err,'-o','Linewidth',1.5);
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