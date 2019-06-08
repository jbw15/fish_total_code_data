
%============calculate the resistive force=================
load('sz2.mat');  % Lateral height of fish 
temp_a = zeros(1,3);
temp_v = zeros(1,3);
temp_body = zeros(1,3);
for nt = 1:num
    for nl = 1:nl_totl
        temp_a(1:2) = accel_s(nt,nl,1:2);
        temp_v(1:2) = vel_int_s(nt,nl,1:2);
        temp_body(1:2) = hh(nt,nl+1,1:2) - hh(nt,nl,1:2);
        temp_a(3) = 0;
        temp_v(3) = 0;
        temp_body(3) = 0;
        accel_unit = temp_a(1:2)/norm(temp_a);
        veloc_unit = temp_v(1:2)/norm(temp_v);
        body_unit = temp_body(1:2)/norm(temp_body);
        cos_theta(nt,nl) = dot(veloc_unit,body_unit)/(norm(veloc_unit)*norm(body_unit));
        if dot(veloc_unit,[-body_unit(2),body_unit(1)])>0
            vy_unit = [-body_unit(2),body_unit(1)];
            
        else 
            vy_unit = [body_unit(2),-body_unit(1)];
       
        end
        acc_amp = norm(temp_a);
        
        f_inert(nt,nl)  = acc_amp*(-dm(nl));               %inertial force
        force_inert(nt,nl,:)=f_inert(nt,nl)*accel_unit;
        
        Cd = 1;dx=0.01;
        Cd2 = 0.06;
        fd1(nt,nl)=-1/2*Cd*dx*sz2(nl)*sum(temp_v.^2)*(1-cos_theta(nt,nl)^2); 
        force_resist2(nt,nl,:) = fd1(nt,nl)*vy_unit; % - fd2(nt,nl)*body_unit;      
    end
    
end
% lateral force

Fy_resist2 =  force_resist2(:,:,2);
Fy_cfd     =  sts_int_sl(:,:,2);
Fx_cfd     =  sts_int_sl(:,:,1);

%  density of lateral force
dx=0.01;
for nt=1:num  
    fy_rea(nt,:)=force_reactive(nt,:)./dx;    %EBT force
    fy_res(nt,:)=Fy_resist2(nt,:)./dx;        %drag force
    fy_cfd(nt,:)=Fy_cfd(nt,:)./dx;
end

% fy_res_distri_fine_f = func_smoothing_in_2D(fy_res, num, nl_totl, refine1, refine2);
% fy_rea_distri_fine_f = func_smoothing_in_2D(fy_rea, num, nl_totl, refine1, refine2);
% fy_cfd_distri_fine_f = func_smoothing_in_2D(fy_cfd, num, nl_totl, refine1, refine2);
% fy_EBT= fy_res_distri_fine_f + fy_rea_distri_fine_f;


% low_v = -2e-3;
% top_v = -low_v;
% fn  = 'fy ebt';
% figure();
% imagesc(xs,xt,fy_res_distri_fine_f);
% axis xy;
% cb = colorbar;
% xlabel('Head<-Position->Tail','FontName','Times');
% ylabel('Time','FontName','Times');
% title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
% set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 0.5 1.0 1.5 2.0]);
% colormap('jet');
% set(cb,'xtick',[low_v 0 top_v])
% caxis([low_v top_v])
% xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
% ylabel('\fontsize{20}\fontname{Times new roman}Time')
% hold off

