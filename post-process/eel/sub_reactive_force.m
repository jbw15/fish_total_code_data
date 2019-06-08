load('dm.mat'); % fish mass distribution
load('dm3.mat'); % added mass distribution



for nt = 1:num
    for nl = 1:nl_totl
        dx = hh(nt,nl+1,1) - hh(nt,nl,1);
        dh_dx(nt,nl) = (hh(nt,nl+1,2)-hh(nt,nl,2))/dx;
        dhdt_dx(nt,nl) = (vel_int_h(nt,nl+1,2)-vel_int_h(nt,nl,2))/dx;
    end    
    for nl = 2:nl_totl-1
        dx_2 = hhs(nt,nl+1,1) - hhs(nt,nl-1,1);
        dhdx_1 = (hhs(nt,nl+1,2)-hhs(nt,nl,2))/(hhs(nt,nl+1,1)-hhs(nt,nl,1));
        dhdx_2 = (hhs(nt,nl,2)-hhs(nt,nl-1,2))/(hhs(nt,nl,1)-hhs(nt,nl-1,1));
        d2hdx2(nt,nl) = (dhdx_1-dhdx_2)/(hh(nt,nl+1,1)-hh(nt,nl,1));
        dmdx(nt,nl) = (dm3(nl+1)-dm3(nl-1))/dx_2;
    end
    nl = 1;
    d2hdx2(nt,nl) = 2*d2hdx2(nt,nl+1)-d2hdx2(nt,nl+2);
    dmdx(nt,nl) = 2*dmdx(nt,nl+1)-dmdx(nt,nl+2);
    
    nl = nl_totl;
    d2hdx2(nt,nl) = 2*d2hdx2(nt,nl-1)-d2hdx2(nt,nl-2);
    dmdx(nt,nl) = 2*dmdx(nt,nl-1)-dmdx(nt,nl-2);
end

for nt = 2:num-1
    accel_s(nt,:,:) = (vel_int_s(nt+1,:,:)-vel_int_s(nt-1,:,:))/(2*dt*n_interval);
end
nt = 1;
accel_s(nt,:,:)  = 2*accel_s(nt+1,:,:) - accel_s(nt+2,:,:);
nt = num;
accel_s(nt,:,:)  = 2*accel_s(nt-1,:,:) - accel_s(nt-2,:,:);

temp = accel_s(:,:,2)+2*U*dhdt_dx+U^2*d2hdx2;

force_reactive1 = -dm3.*temp;
force_reactive2 = -U*dmdx.*(vel_int_s(:,:,2)+U*dh_dx);
force_reactive = force_reactive1 + force_reactive2;

% force_reactive_fine_f = func_smoothing_in_2D(force_reactive, num, nl_totl,refine1, refine2);

% for nt=1:num  
%     fy_rea(nt,:)=force_reactive(nt,:)./ds;
% end
% fy_rea_distri_fine_f = func_smoothing_in_2D(fy_rea, num, nl_totl, refine1, refine2);
% rea_amp = sqrt(2*var(fy_rea_distri_fine_f));

% figure();
% imagesc(xs,xt*2, fy_rea_distri_fine_f);
% axis xy;
% colorbar;
% xlabel('Head<-Position->Tail','FontName','Times');
% ylabel('Time','FontName','Times');
% hold on
% title('reactive force')
% set(gca,'FontSize',20);
% colormap('jet');
% % caxis([-2e-2 2e-2])
% xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
% ylabel('\fontsize{20}\fontname{Times new roman}Time')





        

