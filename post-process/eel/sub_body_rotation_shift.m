
% Computing the rotating angle of eel
for nt = 1:num
    for nl=2:nl_body
        pnt_body1(nt,nl,1:3) = mean(pnt_body(nt,(nl-2)*40+2:(nl-1)*40+1,1:3));
    end
    pnt_body1(nt,1,1:3) = pnt_body(nt,1,1:3);
    pnt_body1(nt,nl_body+1,1:3) = pnt_body(nt,npoint_body,1:3);
end
theta1 = atan((pnt_body1(num,1,2)-pnt_body1(num/2,1,2))/(pnt_body1(num/2,1,1)-pnt_body1(num,1,1)));
theta2 = atan((pnt_body1(num,end,2)-pnt_body1(num/2,end,2))/(pnt_body1(num/2,end,1)-pnt_body1(num,end,1)));
theta3 = atan((pnt_body1(num/2,1,2)-pnt_body1(1,1,2))/(pnt_body1(1,1,1)-pnt_body1(num/2,1,1)));
theta4 = atan((pnt_body1(num/2,end,2)-pnt_body1(1,end,2))/(pnt_body1(1,end,1)-pnt_body1(num/2,end,1)));
theta11 = 0.5*(theta1+theta2);
theta22 = 0.5*(theta3+theta4);
theta = 0.5*(theta11 + theta22);

%Computing lateral displacement of eel
for nt =1:num
    pnt_body2(nt,:,1) = pnt_body1(nt,:,1)*cos(theta) - pnt_body1(nt,:,2)*sin(theta);
    pnt_body2(nt,:,2) = pnt_body1(nt,:,2)*cos(theta) + pnt_body1(nt,:,1)*sin(theta);
end
y_shift = 0.5*(max(pnt_body2(:,end,2))+min(pnt_body2(:,end,2)));

%==========================================================================
%updating the position, vlocity and force of eel
for nt = 1:num
    pnt_body_t(nt,:,1) = pnt_body(nt,:,1)*cos(theta) - pnt_body(nt,:,2)*sin(theta);
    pnt_body_t(nt,:,2) = pnt_body(nt,:,2)*cos(theta) + pnt_body(nt,:,1)*sin(theta) - y_shift;
    vel_body_t(nt,:,1) = vel_body(nt,:,1)*cos(theta) - vel_body(nt,:,2)*sin(theta);
    vel_body_t(nt,:,2) = vel_body(nt,:,2)*cos(theta) + vel_body(nt,:,1)*sin(theta);
    frc_body_t(nt,:,1) = frc_body(nt,:,1)*cos(theta) - frc_body(nt,:,2)*sin(theta);
    frc_body_t(nt,:,2) = frc_body(nt,:,2)*cos(theta) + frc_body(nt,:,1)*sin(theta);
end

 pnt_body_t(:,:,3) = pnt_body(:,:,3);
 frc_body_t(:,:,3) = frc_body(:,:,3);




 
