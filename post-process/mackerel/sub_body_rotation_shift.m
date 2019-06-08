for nt=1:num
    for i=2:50% X
        PP(nt,i,1:3) = mean(pnt_body(nt,(i-2)*40+2:(i-1)*40+1,1:3));% position
    end
    PP(nt,1  ,1:3)=pnt_body(nt,1,1:3);    
    for i=51:58
        PP(nt,i,1:2)=pnt_tail(nt,(i-49)*15+8,1:2);
    end
    for i=59:65
        PP(nt,i,1:2)=pnt_tail(nt,i+100,1:2);
    end
    PP(nt,51:65,3)=0;
end




theta1=atan((PP(400,1,2)-PP(200,1,2))/(PP(200,1,1)-PP(400,1,1)));
theta2=atan((PP(400,end,2)-PP(200,end,2))/(PP(200,end,1)-PP(400,end,1)));
theta3 = atan((PP(num/2,1,2)-PP(1,1,2))/(PP(1,1,1)-PP(num/2,1,1)));
theta4 = atan((PP(num/2,end,2)-PP(1,end,2))/(PP(1,end,1)-PP(num/2,end,1)));
theta11 = 0.5*(theta1+theta2);
theta22 = 0.5*(theta3+theta4);

theta = 0.5*(theta11 + theta22);

for nt=1:num
    pnt_body_t(nt,:,1) = pnt_body(nt,:,1)*cos(theta) - pnt_body(nt,:,2)*sin(theta);
    pnt_body_t(nt,:,2) = pnt_body(nt,:,2)*cos(theta) + pnt_body(nt,:,1)*sin(theta);
    pnt_tail_t(nt,:,1) = pnt_tail(nt,:,1)*cos(theta) - pnt_tail(nt,:,2)*sin(theta);
    pnt_tail_t(nt,:,2) = pnt_tail(nt,:,2)*cos(theta) + pnt_tail(nt,:,1)*sin(theta);
    
    vel_body_t(nt,:,1) = vel_body(nt,:,1)*cos(theta) - vel_body(nt,:,2)*sin(theta);
    vel_body_t(nt,:,2) = vel_body(nt,:,2)*cos(theta) + vel_body(nt,:,1)*sin(theta);
    vel_tail_t(nt,:,1) = vel_tail(nt,:,1)*cos(theta) - vel_tail(nt,:,2)*sin(theta);
    vel_tail_t(nt,:,2) = vel_tail(nt,:,2)*cos(theta) + vel_tail(nt,:,1)*sin(theta);
    
    frc_body_t(nt,:,1) = frc_body(nt,:,1)*cos(theta) - frc_body(nt,:,2)*sin(theta);
    frc_body_t(nt,:,2) = frc_body(nt,:,2)*cos(theta) + frc_body(nt,:,1)*sin(theta);
    frc_tail_t(nt,:,1) = frc_tail(nt,:,1)*cos(theta) - frc_tail(nt,:,2)*sin(theta);
    frc_tail_t(nt,:,2) = frc_tail(nt,:,2)*cos(theta) + frc_tail(nt,:,1)*sin(theta);
    
    sts_tail_t(nt,:,1) = sts_tail(nt,:,1)*cos(theta) - sts_tail(nt,:,2)*sin(theta);
    sts_tail_t(nt,:,2) = sts_tail(nt,:,2)*cos(theta) + sts_tail(nt,:,1)*sin(theta);
end

    pnt_body_t(:,:,3) = pnt_body(:,:,3);
    pnt_tail_t(:,:,3) = pnt_tail(:,:,3);
    frc_body_t(:,:,3) = frc_body(:,:,3);
    frc_tail_t(:,:,3) = frc_tail(:,:,3);
    sts_tail_t(:,:,3) = sts_tail(:,:,3);
    
    Yshift = ( max(pnt_tail_t(:,npoint_tail,2))  + min(pnt_tail_t(:,npoint_tail,2)) )/2;
    pnt_body_t(:,:,2) = pnt_body_t(:,:,2) - Yshift;
    pnt_tail_t(:,:,2) = pnt_tail_t(:,:,2) - Yshift;
    
