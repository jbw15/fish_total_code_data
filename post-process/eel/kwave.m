%calculate for wave volecity
% clear all
% clc

% load('tt.mat')
tt = (1:400)/200;
s=0:0.01:1;
tt = repmat(tt,size(s,2),1);
tt = tt';
s = repmat(s,size(tt,1),1);
tmag=sqrt(2)*std(Torque_totl);
tmag=repmat(tmag,size(tt,1),1);
b0 = [0,0];
[aa,result] = fminsearch(@(xx)sum(sum(abs(tmag.*sin(xx(1)*s-2*pi*tt+xx(2))-Torque_totl))),b0);
aa(1);
vt = 2*pi/aa(1)