function [ ff_fine_f ] = func_adjust_num_segment(ff,dx,ss,mm)
xt0 = (1:400)/200;
xt = (1:500)/250;
numt = length(xt);
if mm == 1
   xs = dx/2:dx:1-dx/2;
else 
   xs = 0:dx:1;
end
% [ss_1,xt_1] = meshgrid(ss,xt0);
[ss_2,xt_2] = meshgrid(xs,xt);

ff_fine_f = interp2(ss,xt0,ff,ss_2,xt_2,'spline');

end

