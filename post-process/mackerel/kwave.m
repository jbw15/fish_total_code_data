%calculate for wave volecity
% clear all
% clc
tt = (1:400)/200;
tt = repmat(tt,65,1);
tt = tt';
load('ss.mat')
s = repmat(ss,400,1);
tmag=sqrt(2)*std(Torque_totl);
tmag=repmat(tmag,size(tt,1),1);
 b0 = [0,0];
[aa,result] = fminsearch(@(xx)sum(sum(abs(tmag.*sin(xx(1)*s-2*pi*tt+xx(2))-Torque_totl))),b0);
aa(1);
a = 2*pi/aa(1)
