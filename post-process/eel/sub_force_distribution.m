ds = 0.01;
for nt = 1:num
    for nl = 2:nl_body-1
        nslt = (nl-2)*nbody_peri+2:nl*nbody_peri+1;
        sts_int_s(nt,nl,1:3) = 0.5*sum(frc_body_t(nt,nslt,1:3))/ds;
        sts_int_sl(nt,nl,1:3) = 0.5*sum(frc_body_t(nt,nslt,1:3));
    end
    nl = 1;
    sts_int_s(nt,1,1:3) = (frc_body_t(nt,1,1:3)+0.5*sum(frc_body_t(nt,2:nbody_peri+1,1:3)))/ds;
    sts_int_sl(nt,1,1:3) = frc_body_t(nt,1,1:3)+0.5*sum(frc_body_t(nt,2:nbody_peri+1,1:3));
    nl = nl_body;
    sts_int_s(nt,nl_body,1:3) = (frc_body_t(nt,npoint_body,1:3)+0.5*sum(frc_body_t(nt,npoint_body-nbody_peri:npoint_body-1,1:3)))/ds;
    sts_int_sl(nt,nl_body,1:3) = frc_body_t(nt,npoint_body,1:3)+0.5*sum(frc_body_t(nt,npoint_body-nbody_peri:npoint_body-1,1:3));
end

%smoothing the sts_int_s
ss = 0:0.01:1;
ms = (ss(1:end-1)+ss(2:end))/2;
for j = 1:2
    sts_temp = sts_int_s(:,:,j);
    sts_temp_f = func_smoothing_in_2D(sts_temp, num, nl_totl,refine1, refine2,ms);
    sts_int_fine_f(:,:,j) = sts_temp_f;
end

index = 2;
if index ==1
    low_v = -4e-3;
    top_v =  4e-3;
    fn = 'Fx';
elseif index ==2
    low_v = -5e-3;
    top_v =  5e-3;
    fn = 'Fy';
end
%=======================================fy=================================
figure();
mass = 0:0.01:1;
xxs = 0.5*(mass(1:end-1)+mass(2:end));
xxt = (1:400)/200;
imagesc(xxs,xxt,sts_int_s(:,:,index));
hold on
xs = linspace(0,1,2*nl_totl); xt = (1:800)/400;
contour(xs(1:end-1),xt,sts_int_fine_f(:,1:end-1,index),[0 0],'--k','LineWidth',3)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
ylabel('Time','FontName','Times','FontSize',20);
% title([fn  ' distribution from CFD'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'FontName','Times','xtick',[0.0 0.5 1.0], 'ytick',[0.0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
xlabel('\fontsize{20}\fontname{Times}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times}Time')
set(gca,'Fontsize',20,'FontName','Times')
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
% gg = gcf;
% print(gg,'result-figure/e-fy-1.eps','-r600','-depsc')
hold off

%=================================fx=======================================
index = 1;
if index ==1
    low_v = -4e-3;
    top_v =  4e-3;
    fn = 'Fx';
elseif index ==2
    low_v = -5e-3;
    top_v =  5e-3;
    fn = 'Fy';
end
figure();
mass = 0:0.01:1;
xxs = 0.5*(mass(1:end-1)+mass(2:end));
xxt = (1:400)/200;
imagesc(xxs,xxt,sts_int_s(:,:,index));
% hold on
% xs = linspace(0,1,2*nl_totl); xt = (1:800)/400;
% contour(xs(1:end-1),xt,sts_int_fine_f(:,1:end-1,index),[1e-10,1e-10],'--k','LineWidth',3)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
ylabel('Time','FontName','Times','FontSize',20);
% title([fn  ' distribution from CFD'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'FontName','Times','xtick',[0.0 0.5 1.0], 'ytick',[0.0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
xlabel('\fontsize{20}\fontname{Times}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times}Time')
set(gca,'Fontsize',20,'FontName','Times')
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
% gg = gcf;
% print(gg,'result-figure/e-fx-1.eps','-r600','-depsc')
hold off

