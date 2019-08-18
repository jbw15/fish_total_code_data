
% Torque calculation integrating from head to tail
for nt = 1:num
    for nl = 2:nl_body
        rx_vector = pnt_body_t(nt,1:(nl-2)*nbody_peri+1,1) - hh(nt,nl,1);
        ry_vector = pnt_body_t(nt,1:(nl-2)*nbody_peri+1,2) - hh(nt,nl,2);
        temp1     = frc_body_t(nt,1:(nl-2)*nbody_peri+1,2).*rx_vector;
        temp2     = frc_body_t(nt,1:(nl-2)*nbody_peri+1,1).*ry_vector;
        torque_fluid1(nt,nl) = sum(temp1) - sum(temp2);
    end
    
    for nl = 1:nl_tail-1
        rx_vector    =  pnt_body_t(nt,1:npoint_body,1) - hh(nt,nl+nl_body,1);
        ry_vector    =  pnt_body_t(nt,1:npoint_body,2) - hh(nt,nl+nl_body,2);
        torque_temp1 =  frc_body_t(nt,1:npoint_body, 2).*rx_vector;
        torque_temp2 =  frc_body_t(nt,1:npoint_body, 1).*ry_vector;
        torque_body  =  sum(torque_temp1) - sum(torque_temp2);
       
        torque_tail = 0;    
        for nn = 2*ntail_vtcl+1:npoint_tail
            if pnt_tail_t(nt,nn,1)<hh(nt,nl+nl_body,1)
                rx_vector    =  pnt_tail_t(nt,nn,1) - hh(nt,nl+nl_body,1);
                ry_vector    =  pnt_tail_t(nt,nn,2) - hh(nt,nl+nl_body,2);
                tq_temp = frc_tail_t(nt,nn, 2)*rx_vector - frc_tail_t(nt,nn, 1)*ry_vector;
                torque_tail  =  torque_tail + tq_temp;
            end
        end
        torque_fluid1(nt,nl+nl_body) = torque_body + torque_tail;
    end
    
    for nl = 1:nl_totl-1
        rx_vector    = hhs(nt,1:nl, 1) - hh(nt,nl,1);
        ry_vector    = hhs(nt,1:nl, 2) - hh(nt,nl,2);
        torque_temp1 = force_inert(nt,1:nl, 2).*rx_vector;
        torque_temp2 = force_inert(nt,1:nl, 1).*ry_vector;
        torque_inert1(nt,nl) = sum(torque_temp1) - sum(torque_temp2);
    end
end
Torque_totl1 = zeros(num,nl_totl);
Torque_totl1(:,2:nl_totl) = torque_fluid1 + torque_inert1;

% Torque calculation integrating from tail to head
for nt = 1:num
    for nl = nl_totl-1:-1:nl_body+1
        torque_tail = 0;
        for nn = 2*ntail_vtcl+1:npoint_tail
            if pnt_tail_t(nt,nn,1)>hh(nt,nl,1)
                rx_vector    =  pnt_tail_t(nt,nn,1) - hh(nt,nl,1);
                ry_vector    =  pnt_tail_t(nt,nn,2) - hh(nt,nl,2);
                tq_temp      = frc_tail_t(nt,nn, 2)*rx_vector - frc_tail_t(nt,nn, 1)*ry_vector;
                torque_tail  =  torque_tail + tq_temp;
            end
        end
        torque_fluid2(nt,nl) = torque_tail;
    end
    
    for nl = nl_body:-1:2
        rx_vector    =  pnt_tail_t(nt,2*ntail_vtcl+1:npoint_tail, 1) - hh(nt,nl,1);
        ry_vector    =  pnt_tail_t(nt,2*ntail_vtcl+1:npoint_tail, 2) - hh(nt,nl,2);
        torque_temp1 =  frc_tail_t(nt,2*ntail_vtcl+1:npoint_tail, 2).*rx_vector;
        torque_temp2 =  frc_tail_t(nt,2*ntail_vtcl+1:npoint_tail, 1).*ry_vector;
        torque_tail  =  sum(torque_temp1) - sum(torque_temp2);
        
        rx_vector    = pnt_body_t(nt,(nl-2)*nbody_peri+1:npoint_body, 1) - hh(nt,nl,1);
        ry_vector    = pnt_body_t(nt,(nl-2)*nbody_peri+1:npoint_body, 2) - hh(nt,nl,2);
        temp1        = frc_body_t(nt,(nl-2)*nbody_peri+1:npoint_body, 2).*rx_vector;
        temp2        = frc_body_t(nt,(nl-2)*nbody_peri+1:npoint_body, 1).*ry_vector;
        torque_body  = sum(temp1) - sum(temp2);
        
        torque_fluid2(nt,nl) = torque_body + torque_tail;
    end
    
    for nl = nl_totl-1:-1:1
        rx_vector    = hhs(nt,nl:nl_totl-1, 1) - hh(nt,nl,1);
        ry_vector    = hhs(nt,nl:nl_totl-1, 2) - hh(nt,nl,2);
        torque_temp1 = force_inert(nt,nl:nl_totl-1, 2).*rx_vector;
        torque_temp2 = force_inert(nt,nl:nl_totl-1, 1).*ry_vector;
        torque_inert2(nt,nl) = sum(torque_temp1) - sum(torque_temp2);
    end
    
end
Torque_totl2 = zeros(num,nl_totl);   
Torque_totl2(:,1:nl_totl-1) = -(torque_fluid2 + torque_inert2);

load('ss.mat');
Torque_totl = (1-ss).*Torque_totl1 + ss.*Torque_totl2;
Torque_totl_fine_f = func_smoothing_in_2D(Torque_totl, num, nl_totl, refine1, refine2,ss);
 



low_v = -4e-4;
top_v = -low_v;
fn  = 'torque';
load('ss.mat')
ss1 = 0:0.01:1;
xxt = (1:400)/200;
% [ss_1,xxt_1] = meshgrid(ss,xxt);
[ss_2,xxt_2] = meshgrid(ss1,xxt);
Torque_totl0 = interp2(ss,xxt,Torque_totl,ss_2,xxt_2,'spline');
figure();
imagesc(ss1,xxt, Torque_totl0);
hold on
plot([ss(1) 0.2205],[1.75 2],'k-','linewidth',2)
plot([ss(1) 0.8152],[1.25 2],'k--' ,'linewidth',2)
plot([ss(1) 1],[0.75 1.75],'k-','linewidth',2)
plot([ss(1) 1],[0.25 1.25],'k--' ,'linewidth',2)
plot([0.2205 1],[0 0.75],'k-','linewidth',2)
plot([0.8152 1],[0 0.25],'k--' ,'linewidth',2)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times','Fontsize',20);
ylabel('Time','FontName','Times','Fontsize',20);
%title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 12)
set(gca,'FontName','Times','FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 0.5 1.0 1.5 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
xlabel('\fontsize{20}\fontname{Times}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times}Time')
gg = gcf;
print(gg,'result-figure/m-torque-1.eps','-r600','-depsc')
hold off

% Torque by viscosity force
scale = 1;
load('ss.mat');
tt = (1:400)/200;
a0 = 1*scale; a1= -3.2*scale; a2 = 5.6*scale;
k = 2*pi/1.0; omega = 2*pi;
A0 = a0 + a1*ss + a2*ss.^2;

for nt = 1:num
        kappa(nt,:)    = A0.*sin(k*ss - omega*tt(nt));
        kappa_dot(nt,:)= -omega*A0.*cos(k*ss - omega*tt(nt));
        kappa_dot(nt,:) = kappa_dot(nt,:);
end
T1 = std(Torque_totl);
dkappa0 = std(kappa_dot);
kk1=0.4*T1./dkappa0;
for nl=1:nl_totl
    torque_vis(:,nl) = kk1(:,nl)*kappa_dot(:,nl);
end
Torque_sum_vis = Torque_totl + torque_vis;
Torque_sum_vis_fine_fs = func_smoothing_in_2D(Torque_sum_vis, num, nl_totl, refine1, refine2,ss);
low_v = -4e-4;
top_v = -low_v;
fn = 'torque by viscosity force';
load('ss.mat')
ss1 = 0:0.01:1;
xxt = (1:400)/200;
% [ss_1,xxt_1] = meshgrid(ss,xxt);
[ss_2,xxt_2] = meshgrid(ss1,xxt);
Torque_sum_vis0 = interp2(ss,xxt,Torque_sum_vis,ss_2,xxt_2,'spline');
figure;
imagesc(ss1,xxt,Torque_sum_vis0);
hold on
plot([ss(1) 0.2205],[1.75 2],'k-','linewidth',2)
plot([ss(1) 0.8152],[1.25 2],'k--' ,'linewidth',2)
plot([ss(1) 1],[0.75 1.75],'k-','linewidth',2)
plot([ss(1) 1],[0.25 1.25],'k--' ,'linewidth',2)
plot([0.2205 1],[0 0.75],'k-','linewidth',2)
plot([0.8152 1],[0 0.25],'k--' ,'linewidth',2)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times','Fontsize',20);
ylabel('Time','FontName','Times','Fontsize',20);
%title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'Fontname','Times','FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
xlabel('\fontsize{20}\fontname{Times}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times}Time')
gg = gcf;
print(gg,'result-figure/m-torque-vis-1.eps','-r600','-depsc')
hold off
% 
