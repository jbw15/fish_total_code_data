
%==========================================================================
%Torque calculation integrating from head to tail
for nt = 1:num
    
%torque by lateral force
    for nl = 2:nl_totl+1
        rx_vector = pnt_body_t(nt,1:(nl-2)*nbody_peri+1,1) - hh(nt,nl,1);
        ry_vector = pnt_body_t(nt,1:(nl-2)*nbody_peri+1,2) - hh(nt,nl,2);
        
        temp1     = frc_body_t(nt,1:(nl-2)*nbody_peri+1,2).*rx_vector;
        temp2     = frc_body_t(nt,1:(nl-2)*nbody_peri+1,1).*ry_vector;
        torque_fluid1(nt,nl) = sum(temp1) - sum(temp2);
    end

%torque by inertial force
    for nl = 2:nl_totl+1
        rx_vector = hhs(nt,1:nl-1,1) - hh(nt,nl,1);
        ry_vector = hhs(nt,1:nl-1,2) - hh(nt,nl,2);
        temp1 = force_inert(nt,1:nl-1,2).*rx_vector;
        temp2 = force_inert(nt,1:nl-1,1).*ry_vector;
        torque_inert1(nt,nl) = sum(temp1) - sum(temp2);
    end
end
Torque_totl1 = torque_fluid1 + torque_inert1;
Torque_totl1(:,1)=0;
%==========================================================================
%Torque calculation integrating from tail to head
for nt = 1:num

%torque by lateral sforce
    rx_vector = pnt_body_t(nt,1:npoint_body,1) - hh(nt,1,1);
    ry_vector = pnt_body_t(nt,1:npoint_body,2) - hh(nt,1,2);
    temp1     = frc_body_t(nt,1:npoint_body,2).*rx_vector;
    temp2     = frc_body_t(nt,1:npoint_body,1).*ry_vector;
    torque_fluid2(nt,1) = sum(temp1) - sum(temp2);
    
    for nl = 2:nl_totl
        rx_vector = pnt_body_t(nt,(nl-2)*nbody_peri+2:npoint_body,1) - hh(nt,nl,1);
        ry_vector = pnt_body_t(nt,(nl-2)*nbody_peri+2:npoint_body,2) - hh(nt,nl,2);
        temp1     = frc_body_t(nt,(nl-2)*nbody_peri+2:npoint_body,2).*rx_vector;
        temp2     = frc_body_t(nt,(nl-2)*nbody_peri+2:npoint_body,1).*ry_vector;
        torque_fluid2(nt,nl) = sum(temp1) - sum(temp2);
    end
    
%torque by inertial force
    for nl = 1:nl_totl
        rx_vector = hhs(nt,nl:nl_totl,1) - hh(nt,nl,1);
        ry_vector = hhs(nt,nl:nl_totl,2) - hh(nt,nl,2);
        temp1 = force_inert(nt,nl:nl_totl,2).*rx_vector;
        temp2 = force_inert(nt,nl:nl_totl,1).*ry_vector;
        torque_inert2(nt,nl) = sum(temp1) - sum(temp2);
    end
end
Torque_totl2 = -(torque_fluid2 + torque_inert2);

Torque_totl2(:,101)=0;
s = 0:0.01:1;
Torque_totl = (1-s).*Torque_totl1+s.*Torque_totl2;
tt = (1:400)/200;
ss = 0:0.01:1;
a_max = 11.41;
k = 2*pi/0.59; omega = 2*pi;
a = a_max.*exp(ss-1);

for nt = 1:num
    kappa(nt,:) = a.*sin(k*ss-omega*tt(nt));
end

low_v = -2e-4;
top_v = -low_v;
fn  = 'torque';
figure();
xxt = (1:400)/200;ss = 0:0.01:1;
imagesc(ss,xxt,Torque_totl);
hold on
plot([0.75 1],[0 0.4450],'k-' ,'linewidth',2)
plot([0.45 1],[0 0.9450],'k--' ,'linewidth',2)
plot([0.16 1],[0 1.4450],'k-' ,'linewidth',2)
plot([0 1],[0.25 1.9450],'k--' ,'linewidth',2)
plot([0 0.75],[0.75 2],'k-' ,'linewidth',2)
plot([0 0.45],[1.25 2],'k--' ,'linewidth',2)
plot([0 0.16],[1.75 2],'k-' ,'linewidth',2)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 0.5 1.0 1.5 2.0]);
set(gca,'yticklabel',{'0','0.5','1.0','1.5','2.0'});
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
hold off

tt = (1:400)/200;
ss = 0:0.01:1;
a_max = 11.41;
k = 2*pi/0.59; omega = 2*pi;
a = a_max.*exp(ss-1);

for nt = 1:num
    kappa(nt,:) = a.*sin(k*ss-omega*tt(nt));
    kappa_dot(nt,:) = -omega*a.*cos(k*ss-omega*tt(nt));
    kappa_dot(nt,:) = kappa_dot(nt,:);
end        

T1 = std(Torque_totl);
dkappa0 = std(kappa_dot);
kk1=0.4*T1./dkappa0;
for nl=1:nl_totl+1
    torque_vis(:,nl) = kk1(:,nl)*kappa_dot(:,nl);
end
Torque_sum_vis = Torque_totl + torque_vis;
low_v = -2e-4;
top_v = -low_v;
xxt = (1:400)/200;ss = 0:0.01:1;
figure;
imagesc(ss,xxt,Torque_sum_vis);
hold on
plot([0.75 1],[0 0.4450],'k-' ,'linewidth',2)
plot([0.45 1],[0 0.9450],'k--' ,'linewidth',2)
plot([0.16 1],[0 1.4450],'k-' ,'linewidth',2)
plot([0 1],[0.25 1.9450],'k--' ,'linewidth',2)
plot([0 0.75],[0.75 2],'k-' ,'linewidth',2)
plot([0 0.45],[1.25 2],'k--' ,'linewidth',2)
plot([0 0.16],[1.75 2],'k-' ,'linewidth',2)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
hold off
       
        
        
        
        
