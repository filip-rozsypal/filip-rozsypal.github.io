close all; clear all;clc;

time = 0.5:0.5:75;

A_0 = 1;
L_0 = 3;
K0 = 90;

g = 0.03;
n = 0.05;

A = A_0*exp(g*time);
L = L_0*exp(n*time);

figure(1)
plot(time,[A;L]','linewidth',2);
legend('A(t)','L(t)',2);legend('boxoff')
xlabel('time');






figure(2)
plot(time,[log(A);log(L)]','linewidth',2)
legend('log(A(t))','log(L(t))',2);legend('boxoff');
xlabel('time');





capital = 0:0.05:100;
alpha = 0.6;
delta = 0.1;
s = 0.7;

f = @(x) x.^alpha;

figure(3)
plot(capital,[s*f(capital);(delta+n+g)*capital]','linewidth',2)
xlabel('k'); legend('sf(k)','(\delta+n+g)k',2);legend('boxoff');
axis([0 40 0 8])

figuresize(7,6,'centimeters')
target=strcat(target_dir,'kdot','.pdf');
print(gcf, '-dpdf', '-r300', target);


figure(4)
plot(capital,[s*f(capital)./capital;(delta+n+g)*ones(size(capital))]','linewidth',2)
xlabel('k'); %ylabel('$\dot{k}/k$','interpreter','latex')
legend('sf(k)/k','\delta+n+g',1);legend('boxoff');
axis([1 40 0 0.7])




figure(5)
plot(capital,s*f(capital)./capital-(delta+n+g)*ones(size(capital)),'linewidth',2)
xlabel('k'); ylabel('$\dot{k}/k$','interpreter','latex')
legend('sf(k)/k-(\delta+n+g)',1);legend('boxoff');
axis([0 40 -0.1 0.5])
hline(0,'-k');




figure(6)
plot(capital,s*f(capital)-(delta+n+g)*capital,'linewidth',2)
xlabel('k'); ylabel('$\dot{k}}$','interpreter','latex')
legend('sf(k)-(\delta+n+g)k',1);legend('boxoff');
axis([0 40 -0.5 2])
hline(0,'-k');






figure(7)
plot(capital,s*f(capital)-(delta+n+g-1)*capital,'linewidth',2);hold on;
plot(capital,capital,'--k','linewidth',1.5);hold off;
xlabel('k_t'); ylabel('$k_{t+1}$','interpreter','latex')
legend('sf(k)-(\delta+n+g-1)k','45 degree line',2);legend('boxoff');
axis([0 35 0 35])
hline(0,'-k');




%%

dk_star=@(x) s*f(x)-(delta+n+g)*x;
k_star=fsolve(dk_star,100); %get the steady state capital

%% convergence

end_period=100;
time = 1:end_period;

A = zeros(size(time));
L = A;
K = L;
k = K;

A(1) = 1;
L(1) = 3;
K(1) =K0;
k(1) = K(1)/(A(1)*L(1));



for t=2:end_period
    A(t) = (1+g)*A(t-1);
    L(t) = (1+n)*L(t-1);
    k(t) = s*f(k(t-1))-(delta+n+g-1)*k(t-1);
end

K = A.*L.*k;
y = k.^alpha;
Y = A.*L.*y;
c = (1-s)*y;
C = (1-s)*Y;


%%
figure(9)
plot(time,[k],'linewidth',2);
legend('k',2); legend('boxoff');
hline(k_star);
xlabel('time')

figuresize(7,6,'centimeters')
target=strcat(target_dir,'k_time','.pdf');
print(gcf, '-dpdf', '-r300', target);


figure(8)
plot(time(2:end),(k(2:end)-k(1:end-1))./k(1:end-1),'linewidth',2);
legend('growth rate of k',2); legend('boxoff');
xlabel('time')

figuresize(7,6,'centimeters')
target=strcat(target_dir,'k_gr_time','.pdf');
print(gcf, '-dpdf', '-r300', target);

%%



figure(10)
plot(time,[k;y;c]);
legend('k','y','c',2); legend('boxoff');
hline(k_star);
xlabel('time')

figure(11)
plot(time,[y;c]);
legend('y','c',2); legend('boxoff');
hline(k_star);
xlabel('time')
axis([1 end_period 0 2.5])


%%
figure(12)
plot(time,[Y;C])
axis([1 end_period 0 10000])
xlabel('time');
legend('Y','C',2); legend('boxoff');

figure(13)
plot(time,[log(Y);log(C)])
axis([1 end_period/5 0 5])
xlabel('time');
legend('log(Y)','log(C)',2); legend('boxoff');


figure(14)
plot(time,[log(Y);log(Y./L);log(Y./(A.*L))],'linewidth',2)
axis([0 end_period 1.5 8])
xlabel('time'); 
legend('log(Y)','log(Y/L)','log(Y/AL)',2); legend('boxoff');

figuresize(7,6,'centimeters')
target=strcat(target_dir,'Y_gr_time','.pdf');
print(gcf, '-dpdf', '-r300', target);


%% Change in n
n_old = n;
n_new = n/5;


end_period=100;
change_period = 50;

time = 1:end_period;

n_time = [n_old*ones(1,change_period-1) n_new*ones(1,end_period-change_period+1)];

%% change in n

A = zeros(size(time));
L = A;
K = L;
k = K;

A(1) = 1;
L(1) = 3;
K(1) = 20;
k(1) = K(1)/(A(1)*L(1));



for t=2:end_period
    A(t) = (1+g)*A(t-1);
    L(t) = (1+n_time(t))*L(t-1);
    k(t) = s*f(k(t-1))-(delta+n_time(t)+g-1)*k(t-1);
end

K = A.*L.*k;
y = k.^alpha;
Y = A.*L.*y;
c = (1-s)*y;
C = (1-s)*Y;

%%

figure(20)
plot(1941:2040,log(L),'linewidth',2);
vline(1989);
legend('log(N)',2);legend('boxoff')
set(gca,'XTick',1950:20:2050);
axis([1950 2040 1 5])




figure(21)
hl=plot(capital,[s*f(capital);(delta+n_old+g).*capital;(delta+n_new+g).*capital]','linewidth',2);
xlabel('k'); %ylabel('$\dot{k}/k$','interpreter','latex')
legend('sf(k)','(\delta+n_{old}+g)k','(\delta+n_{new}+g)k',4);legend('boxoff');
set(hl(2),'linewidth',1.5,'color','k','lineStyle','--');
axis([0 70 0 10])





figure(22)
hl=plot(1941:2040,k,'linewidth',2);
xlabel('time'); %ylabel('$\dot{k}/k$','interpreter','latex')
legend('k',2);legend('boxoff');
set(gca,'XTick',1950:20:2050);
axis([1940 2040 0 70])
vline(1989);




figure(23)
plot(1941:2040,[log(Y);log(Y./L);log(Y./(A.*L))],'linewidth',2)
set(gca,'XTick',1950:20:2050);
axis([1940 2040 1 8])
xlabel('time'); vline(1989);
legend('log(Y)','log(Y/L)','log(Y/AL)',2); legend('boxoff');





figure(24)
hl=plot(1942:2040,(k(2:end)-k(1:end-1))./k(1:end-1),'linewidth',2);
xlabel('time'); %ylabel('$\dot{k}/k$','interpreter','latex')
legend('growth rate of k',2);legend('boxoff');
set(gca,'XTick',1950:20:2050);
axis([1940 2040 0 0.15])
vline(1989);










