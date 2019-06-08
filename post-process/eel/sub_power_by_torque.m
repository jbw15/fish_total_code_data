
tt = (1:400)/200;
ss = 0:0.01:1;
a_max = 11.41;
k = 2*pi/0.59; omega = 2*pi;
a = a_max.*exp(ss-1);

for nt = 1:num
    kappa(nt,:) = a.*sin(k*ss-omega*tt(nt));
    kappa_dot(nt,:) = -omega*a.*cos(k*ss-omega*tt(nt));
end

for nt = 1:num
    for nl = 1:nl_totl+1
        power_fluid(nt,nl) = torque_fluid1(nt,nl)*kappa_dot(nt,nl);
        power_inert(nt,nl) = torque_inert1(nt,nl)*kappa_dot(nt,nl);
        power_tt(nt,nl)    = Torque_totl(nt,nl)*kappa_dot(nt,nl);
      
    end
end
power_totl = power_fluid + power_inert;
 
low_v = -6e-3;
top_v = -low_v;
fn    = 'power by torque';
figure;
xxt = (1:400)/200;ss = 0:0.01:1;
imagesc(ss,xxt,power_tt);
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
hold off


power_positive = zeros(nt,nl_totl+1);
for nt = 1:num
    for nl = 1:nl_totl+1
        if power_tt(nt,nl)<0
            power_positive(nt,nl) = 0;
        else
            power_positive(nt,nl) = power_tt(nt,nl);
        end
    end
end
T0 = std(Torque_totl);
kapa0=std(kappa);
% kk = 0.8*T0./kapa0;%2019.4.19查看系数的影响
kk=0.4*T0./kapa0;
for i=1:nl_totl+1
    Torque_e(:,i)=kk(i)*kappa(:,i);
end
Torque_sum=Torque_totl+Torque_e;
for nt = 1:num
    for nl = 1:nl_totl+1
        power_sum(nt,nl)    = Torque_sum(nt,nl)*kappa_dot(nt,nl); 
        power_sum_vis(nt,nl) = Torque_sum_vis(nt,nl)*kappa_dot(nt,nl);
    end
end
power_sum_positive = zeros(nt,nl_totl+1);
for nt = 1:num
    for nl = 1:nl_totl+1
        if power_sum(nt,nl)<0
            power_sum_positive(nt,nl) = 0;
        else
            power_sum_positive(nt,nl) = power_sum(nt,nl);
        end
    end
end
%=====================Torque_sum===========================================
ss = 0:0.01:1;
Torque_totl_fine_fs = func_smoothing_in_2D(Torque_sum, num, nl_totl+1, refine1, refine2,ss);
low_v = -1.5e-4;
top_v = -low_v;
fn    = 'total torque';
xxx = linspace(0,1,(nl_totl+1)*2);
xt = (1:800)/400;
xxt = (1:400)/200;
figure;
imagesc(ss,xxt,Torque_sum);
hold on
contour(xxx(30:165),xt,Torque_totl_fine_fs(:,30:165),[0,0],'--k','LineWidth',3)
% hold on
% plot([0.75 1],[0 0.4450],'k-' ,'linewidth',2)
% plot([0.45 1],[0 0.9450],'k--' ,'linewidth',2)
% plot([0.16 1],[0 1.4450],'k-' ,'linewidth',2)
% plot([0 1],[0.25 1.9450],'k--' ,'linewidth',2)
% plot([0 0.75],[0.75 2],'k-' ,'linewidth',2)
% plot([0 0.45],[1.25 2],'k--' ,'linewidth',2)
% plot([0 0.16],[1.75 2],'k-' ,'linewidth',2)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
hold off     

%===========================power_sum_vis==================================
low_v = -6e-3;
top_v = -low_v;
fn    = 'power_sum_vis';
xxx = linspace(0,1,(nl_totl+1)*2);
xt = (1:800)/400;
xxt = (1:400)/200;
figure;
imagesc(ss,xxt,power_sum_vis);
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
hold off  

%=====================power_sum============================================
low_v = -6e-3;
top_v = -low_v;
fn    = 'power_sum';
xxx = linspace(0,1,(nl_totl+1)*2);
xt = (1:800)/400;
xxt = (1:400)/200;
figure;
imagesc(ss,xxt,power_sum);
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
hold off  

ds=0.01;
for nt = 1:num
    power_h1(nt) = sum(power_tt(nt,:)*ds);
    power_h2(nt) = sum(power_positive(nt,:)*ds);
end
sum(power_h1*dt*10)/2
sum(power_h2*dt*10)/2

% trapz(tt,power_h1)/2
% trapz(tt,power_h2)/2