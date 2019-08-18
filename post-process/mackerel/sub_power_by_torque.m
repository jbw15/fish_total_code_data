
scale = 1;
tt = (1:400)/200;
load('ss.mat');
a0 = 1*scale; a1= -3.2*scale; a2 = 5.6*scale;
k = 2*pi/1.0; omega = 2*pi;
A0 = a0 + a1*ss + a2*ss.^2;

for nt = 1:num
        kappa(nt,:)    = A0.*sin(k*ss - omega*tt(nt));
        kappa_dot(nt,:)= -omega*A0.*cos(k*ss - omega*tt(nt));
end

for nt = 1:num
    for nl = 1:nl_totl
        power_tt(nt,nl)     =  Torque_totl(nt,nl)*kappa_dot(nt,nl);
    end
end

 power_totl_fine_f  = func_smoothing_in_2D(power_tt, num, nl_totl, refine1, refine2,ss);

low_v = -2e-3;
top_v = -low_v;
fn    = 'power';
load('ss.mat')
ss1 = 0:0.01:1;
xxt = (1:400)/200;
% [ss_1,xxt_1] = meshgrid(ss,xxt);
[ss_2,xxt_2] = meshgrid(ss1,xxt);
power_tt0 = interp2(ss,xxt,power_tt,ss_2,xxt_2,'spline');
figure;
imagesc(ss1,xxt,power_tt0);
hold on
xt = (1:800)/400; xs = linspace(0,1,2*nl_totl);
contour(xs(33:113),xt,Torque_totl_fine_f(:,33:113),[0 0],'--k','LineWidth',3)
hold on
plot(ss(29),0.3,'.k','Markersize',30)
hold on
plot(ss(42),0.3,'.k','Markersize',30)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
ylabel('Time','FontName','Times','FontSize',20);
%title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize',12)
set(gca,'FontSize',20,'Fontname','Times','xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
xlabel('\fontsize{20}\fontname{Times}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times}Time')
gg = gcf;
print(gg,'result-figure/m-power-1.eps','-r600','-depsc')
hold off

power_positive = zeros(nt,nl_totl);
for nt = 1:num
    for nl = 1:nl_totl
        if power_tt(nt,nl)<0
            power_positive(nt,nl) = 0;
        else
            power_positive(nt,nl) = power_tt(nt,nl);
        end
    end
end
T0 = std(Torque_totl);
kapa0=std(kappa);
% kk=0.8*T0./kapa0;%查看系数的作用
kk=0.4*T0./kapa0;
for i=1:nl_totl
    Torque_e(:,i)=kk(i)*kappa(:,i);
end
Torque_sum=Torque_totl+Torque_e;
for nt = 1:num
    for nl = 1:nl_totl
        power_sum(nt,nl)    = Torque_sum(nt,nl)*kappa_dot(nt,nl); 
        power_sum_vis(nt,nl) = Torque_sum_vis(nt,nl)*kappa_dot(nt,nl);
    end
end
power_sum_positive = zeros(nt,nl_totl);
for nt = 1:num
    for nl = 1:nl_totl
        if power_sum(nt,nl)<0
            power_sum_positive(nt,nl) = 0;
        else
            power_sum_positive(nt,nl) = power_sum(nt,nl);
        end
    end
end
Torque_totl_fine_fs = func_smoothing_in_2D(Torque_sum, num, nl_totl, refine1, refine2,ss);

%==========================Torque_sum======================================
low_v = -2e-3;
top_v = -low_v;
low_v = -3e-4;
top_v = -low_v;
fn    = 'total torque';
load('ss.mat')
ss1 = 0:0.01:1;
xxt = (1:400)/200;
% [ss_1,xxt_1] = meshgrid(ss,xxt);
[ss_2,xxt_2] = meshgrid(ss1,xxt);
Torque_sum0 = interp2(ss,xxt,Torque_sum,ss_2,xxt_2,'spline');
xs = linspace(0,1,2*nl_totl);xt = (1:800)/400;
figure;
imagesc(ss1,xxt,Torque_sum0);
% hold on
% plot([ss(1) 0.2205],[1.75 2],'k-','linewidth',2)
% plot([ss(1) 0.8152],[1.25 2],'k--' ,'linewidth',2)
% plot([ss(1) 1],[0.75 1.75],'k-','linewidth',2)
% plot([ss(1) 1],[0.25 1.25],'k--' ,'linewidth',2)
% plot([0.2205 1],[0 0.75],'k-','linewidth',2)
% plot([0.8152 1],[0 0.25],'k--' ,'linewidth',2)
hold on
contour(xs(27:101),xt,Torque_totl_fine_fs(:,27:101),[0 0],'--k','LineWidth',3)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
ylabel('Time','FontName','Times','FontSize',20);
%title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'FontName','Times','xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
xlabel('\fontsize{20}\fontname{Times}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times}Time')
gg = gcf;
print(gg,'result-figure/m-torque-ela-1.eps','-r600','-depsc')
hold off

%=============================power_sum====================================
low_v = -2e-3;
top_v = -low_v;
fn    = 'power_sum';
load('ss.mat')
ss1 = 0:0.01:1;
xxt = (1:400)/200;
% [ss_1,xxt_1] = meshgrid(ss,xxt);
[ss_2,xxt_2] = meshgrid(ss1,xxt);
power_sum0 = interp2(ss,xxt,power_sum,ss_2,xxt_2,'spline');
xs = linspace(0,1,2*nl_totl);xt = (1:800)/400;
figure;
imagesc(ss1,xxt,power_sum0);
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
ylabel('Time','FontName','Times','FontSize',20);
% title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'FontName','Times','xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
xlabel('\fontsize{20}\fontname{Times}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times}Time')
gg = gcf;
print(gg,'result-figure/m-power-ela-1.eps','-r600','-depsc')
hold off
% 
% %===========================power_sum_vis==================================
low_v = -2e-3;
top_v = -low_v;
fn    = 'power_sum_vis';
load('ss.mat')
ss1 = 0:0.01:1;
xxt = (1:400)/200;
% [ss_1,xxt_1] = meshgrid(ss,xxt);
[ss_2,xxt_2] = meshgrid(ss1,xxt);
power_sum_vis0 = interp2(ss,xxt,power_sum_vis,ss_2,xxt_2,'spline');
xs = linspace(0,1,2*nl_totl);xt = (1:800)/400;
figure;
imagesc(ss1,xxt,power_sum_vis0);
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
ylabel('Time','FontName','Times','FontSize',20);
% title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'FontName','Times','xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
xlabel('\fontsize{20}\fontname{Times}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times}Time')
gg = gcf;
print(gg,'result-figure/m-power-vis-1.eps','-r600','-depsc')
hold off


load('ss.mat');
dss(2:nl_totl) = ss(2:end)-ss(1:end-1);
dss(1) = 2*dss(2)-dss(3);
for nt = 1:num
    power_h1(nt) = sum(power_tt(nt,:).*dss);
    power_h2(nt) = sum(power_positive(nt,:).*dss);
end
sum(power_h1*10*dt)/2
sum(power_h2*10*dt)/2

% trapz(tt,power_h1)/2
% trapz(tt,power_h2)/2

