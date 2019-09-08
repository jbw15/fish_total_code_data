%Á¦£¬Á¦¾Ø²åÖµ
Torque_125_2 = func_adjust_num_segment(Torque_totl,0.008,ss,2);
fx_125_2 = func_adjust_num_segment(sts_int_s(:,:,1),0.008,ms,1);
fy_125_2 = func_adjust_num_segment(sts_int_s(:,:,2),0.008,ms,1);
save Torque_125_2 Torque_125_2
save fx_125_2 fx_125_2
save fy_125_2 fy_125_2