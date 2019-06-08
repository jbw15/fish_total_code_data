set(0,'defaultfigurecolor','w')
clear all
clc
npoint=14990;
fn =  ['nlist_eel'];
fid = fopen(fn,'r');
nline = 0;
i=0;
while 1
    line = fgetl(fid);
    if(line==-1)
        break;
    end
    str = strfind(line, 'NODE');
    if (~isempty(str))
        for nn=i*20+1:(i+1)*20
            line = fgetl(fid);
            if nn>npoint
                break
            end
            [P(nn,1:4),count1] = sscanf(line, '%f');      % data at each time frame
        end
        i=i+1;
    end
end
pnt(:,1)=P(:,4)/(max(P(:,4))-min(P(:,4)));
pnt(:,2)=P(:,2)/(max(P(:,4))-min(P(:,4)));
pnt(:,3)=P(:,3)/(max(P(:,4))-min(P(:,4)));
load('tri_eel.mat');
dx=0.01;
mx=0.005:dx:0.995;
NX=length(mx);
for nx=1:NX
    my1(nx)=0;my2(nx)=0;
    mz1(nx)=0;mz2(nx)=0;
    for np=1:npoint
        if pnt(np,1)>(mx(nx)-0.5*dx)&&pnt(np,1)<=(mx(nx)+0.5*dx)
            yy=pnt(np,2);
            zz=pnt(np,3);
            if yy>my1(nx)    my1(nx)=yy;     end
            if yy<my2(nx)    my2(nx)=yy;     end
            if zz>mz1(nx)    mz1(nx)=zz;     end
            if zz<mz2(nx)    mz2(nx)=zz;     end
        end
    end
end
% my1 = my1/2; my2 = my2/2; mz1 = mz1/2; mz2 = mz2/2;
% mz1 = mz1/3; mz2 = mz2/3;
sz2=mz1-mz2;
rho=1;
dm=rho*dx*pi/4*(my1-my2).*(mz1-mz2);
%3/22 the first segment and the last are cone
dm(1) = dm(1)/3; dm(100) = dm(100)/3;
dm3=rho*dx*pi/4*sz2.^2;
m=sum(dm);
dt = 5e-4; t = dt:dt:1;num = length(t);
ds = 0.01; s = 0:ds:1; nl_totl = length(s) - 1; ms = ds/2:ds:1-ds/2;
omega = 2*pi; lambda =0.59; k = 2*pi/lambda;
am = 11.41;
a = am*exp(s-1);
nw=5;

for nt = 1:num
    t = nt*dt;
    curv(nt,:) = a.*sin(k*s-omega*t);
    
%     dcurv(nt,:) = -a*omega*sin(k*s+omega*t);
%     oangle(nt,:)=cumsum(curv(nt,:))*ds;
%     x0(nt,:)=cumsum(cos(oangle(nt,:))*ds);
%     y0(nt,:)=cumsum(sin(oangle(nt,:))*ds);
    for nl = 2:nl_totl+1
        
%         xx = [s1,s(1:nl)];the = [ -a*sin(k*s1+omega*t),curv(nt,1:nl)];
         xx = s(1:nl);the = curv(nt,1:nl);
%         dthe = [-a*omega*sin(k*s1+omega*t),dcurv(nt,1:nl)];

        theta(nt,nl) = trapz(xx,the);
%         dtheta(nt,nl) = trapz(xx,dthe);
%         sxx = [0,cos(theta(nt,1:nl))];syy = [0,sin(theta(nt,1:nl))];
         sxx = cos(theta(nt,1:nl));syy = sin(theta(nt,1:nl));
        x0(nt,nl)=trapz(xx,sxx);
        y0(nt,nl)=trapz(xx,syy);
    end
    X0(nt,:)=0.5*(x0(nt,1:nl_totl)+x0(nt,2:nl_totl+1));
    Y0(nt,:)=0.5*(y0(nt,1:nl_totl)+y0(nt,2:nl_totl+1));    
    cx(nt) = sum(dm.*X0(nt,:))/m;
    cy(nt) = sum(dm.*Y0(nt,:))/m;
    
    x1(nt,:)=x0(nt,:)-cx(nt);
    y1(nt,:)=y0(nt,:)-cy(nt);
    X1(nt,:)=X0(nt,:)-cx(nt);
    Y1(nt,:)=Y0(nt,:)-cy(nt);
    cx1(nt)=sum(dm.*X1(nt,:))/m;
    cy1(nt)=sum(dm.*Y1(nt,:))/m;
end
for nt = 1:num
    X1(nt,:)=X1(nt,:)-cx1(nt);
    Y1(nt,:)=Y1(nt,:)-cy1(nt);
    x1(nt,:)=x1(nt,:)-cx1(nt);
    y1(nt,:)=y1(nt,:)-cy1(nt);
    ic1(nt,:)=dm.*(X1(nt,:).^2+Y1(nt,:).^2);
    Icom1(nt)=sum(ic1(nt,:));
    r_2(nt,:) = X1(nt,:).^2 + Y1(nt,:).^2;
    r(nt,:) = sqrt(r_2(nt,:));   
end
for nt = 2:num-1
    vx(nt,:) = (X1(nt+1,:) - X1(nt-1,:))/(2*dt);
    vy(nt,:) = (Y1(nt+1,:) - Y1(nt-1,:))/(2*dt);
    vr(nt,:) = (r(nt+1,:) - r(nt-1,:))/(2*dt);
end
nt = 1;
vx(nt,:) = 2*vx(nt+1,:) - vx(nt+2,:); vy(nt,:) = 2*vy(nt+1,:)-vy(nt+2,:);vr(nt,:) = 2*vr(nt+1,:)-vr(nt+2,:);
nt = num;
vx(nt,:) = 2*vx(nt-1,:) - vx(nt-2,:); vy(nt,:) = 2*vy(nt-1,:)-vy(nt-2,:);vr(nt,:) = 2*vr(nt-1,:)-vr(nt-2,:);
for nt = 1:num
    for nl = 1:nl_totl
        r1(nl,:) = [X1(nt,nl),Y1(nt,nl),0];
        LA(nl) = norm(r1(nl,:));
        if nt<num
           r2(nl,:) = [X1(nt+1,nl),Y1(nt+1,nl),0];
           LB(nl) = norm(r2(nl,:));
        else 
           r2(nl,:) = [X1(1,nl),Y1(1,nl),0];
           LB(nl) = norm(r2(nl,:));
        end
        CS(nl,:)=cross(r1(nl,:),r2(nl,:));
        if CS(nl,3)>0
            sign(nt,nl) = 1;
        else
            sign(nt,nl) = -1;
        end
        DC(nl)=dot(r1(nl,:),r2(nl,:));
        if DC(nl)>LA(nl)*LB(nl)
            DC(nl)=LA(nl)*LB(nl);
        end
    end
    dtheta1(nt,:)=acos(DC./(LA.*LB)).*sign(nt,:);
    w1(nt,:)=dtheta1(nt,:)/dt;
    iw1(nt,:)=ic1(nt,:).*w1(nt,:);
    IW1(nt) = sum(iw1(nt,:));
end
for nt=2:num
    w(nt)=IW1(nt-1)/(Icom1(nt-1));
    %     rot0(nt)=sum(w(3:nt)*dt);
end

for nt = 1:num
    v_2(nt,:) = vx(nt,:).^2 + vy(nt,:).^2;
    v(nt,:) = sqrt(v_2(nt,:));
    vthe(nt,:) = sign(nt,:).*sqrt(v_2(nt,:)-vr(nt,:).^2);  
end
for nt = 2:num-1
    ay(nt,:) = (vy(nt+1,:)-vy(nt-1,:))/(2*dt);
    ax(nt,:) = (vx(nt+1,:)-vx(nt-1,:))/(2*dt);
    athe(nt,:) = (vthe(nt+1,:)-vthe(nt-1,:))/(2*dt);
end
nt = 1;
ay(nt,:) = 2*ay(nt+1,:)-ay(nt+2,:);athe(nt,:) = 2*athe(nt+1,:)-athe(nt+2,:);
ax(nt,:) = 2*ax(nt+1,:)-ax(nt+2,:);
nt = num;
ay(nt,:) = 2*ay(nt-1,:)-ay(nt-2,:);athe(nt,:) = 2*athe(nt-1,:)-athe(nt-2,:);
ax(nt,:) = 2*ax(nt-1,:)-ax(nt-2,:);
for nt = 1:num
    c1(nt,:) = 2.*dm.*(X1(nt,:).*vx(nt,:)+Y1(nt,:).*vy(nt,:));
    c2(nt,:) = dm.*r_2(nt,:);
    c3(nt,:) = dm.*(X1(nt,:).*ay(nt,:)-Y1(nt,:).*ax(nt,:));
end

for nt = 1:num
    C1(nt) = sum(c1(nt,:));
    C2(nt) = sum(c2(nt,:));
    C3(nt) = sum(c3(nt,:));
end
C3 = C3';
for nt = 1:num
    t = nt*dt;
    for n = 1:nw
        cosw(nt,n) = cos(n*omega*t);
        sinw(nt,n) = sin(n*omega*t);
        ncosw(nt,n) = n*omega*cos(n*omega*t);
        nsinw(nt,n) = n*omega*sin(n*omega*t);
    end
end
for nt = 1:num
    t = nt*dt;
    A(nt,1:nw) = C1(nt)*cosw(nt,:) - C2(nt)*nsinw(nt,:);
    A(nt,nw+1:2*nw) = C1(nt)*sinw(nt,:) + C2(nt)*ncosw(nt,:);

end 
 b0 = [-0.2858    0.0000    0.0016   -0.0000  0  2.5334    0.0000   -0.0008    0.0000 0];


options = optimset('MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-10);

[coe,result,~,report] = fminsearch(@(B)myfun(A,B,C3,nw,dt),b0,options);

 for n = 1:num
    wb(n)=0;
    for nn= 1:nw
        wb(n) = wb(n) + coe(nn)*cosw(n,nn);
    end
    for nn=1:nw
        wb(n) = wb(n) + coe(nn+nw)*sinw(n,nn);
    end
end

rot0=cumsum(wb)*dt;
rot=rot0-mean(rot0);drot=diff(rot)/dt;ddrot=diff(drot)/dt;
for nt = 1:num
    %x1(nt,:)=x1(nt,:)-cx1(nt);
    %y1(nt,:)=y1(nt,:)-cy1(nt);    
    x2(nt,:) = x1(nt,:)*cos(rot(nt)) - y1(nt,:)*sin(rot(nt));
    y2(nt,:) = y1(nt,:)*cos(rot(nt)) + x1(nt,:)*sin(rot(nt));   
    X2(nt,:) = X1(nt,:)*cos(rot(nt)) - Y1(nt,:)*sin(rot(nt));
    Y2(nt,:) = Y1(nt,:)*cos(rot(nt)) + X1(nt,:)*sin(rot(nt));    
    cx2(nt)=sum(dm.*X2(nt,:))/m;
    cy2(nt)=sum(dm.*Y2(nt,:))/m;
end
for nt=1:num
    X2(nt,:)=X2(nt,:)-cx2(nt);
    Y2(nt,:)=Y2(nt,:)-cy2(nt);
    
    x2(nt,:)=x2(nt,:)-cx2(nt);
    y2(nt,:)=y2(nt,:)-cy2(nt);
    
    ic2(nt,:)=dm.*(X2(nt,:).^2+Y2(nt,:).^2);
    Icom2(nt)=sum(ic2(nt,:));
end
aicom = mean(Icom2);
for i=2:num-1
    DI(i)=(Icom2(i+1)-Icom2(i-1))/(2*dt);
end
DI(1)=2*DI(2)-DI(3);
DI(num)=2*DI(num-1)-DI(num-2);
for nt = 2:num-1
       vel_int_s(nt,:,1) = (X2(nt+1,:)-X2(nt-1,:))/(2*dt);  
       vel_int_s(nt,:,2) = (Y2(nt+1,:)-Y2(nt-1,:))/(2*dt);
end
vel_int_s(1,:,:) = 2*vel_int_s(2,:,:) - vel_int_s(3,:,:);
vel_int_s(num,:,:) = 2*vel_int_s(num-1,:,:) - vel_int_s(num-2,:,:);

for nt = 2:num-1
    accel_s(nt,:,:) = (vel_int_s(nt+1,:,:)-vel_int_s(nt-1,:,:))/(2*dt);
end
nt = 1;
accel_s(nt,:,:)  = 2*accel_s(nt+1,:,:) - accel_s(nt+2,:,:);
nt = num;
accel_s(nt,:,:)  = 2*accel_s(nt-1,:,:) - accel_s(nt-2,:,:);

for nt = 1:num

        for nl = 2:nl_totl+1
        rx_vector = X2(nt,1:nl-1) - x2(nt,nl); 
        ry_vector = Y2(nt,1:nl-1) - y2(nt,nl);
        temp1 = -dm(1:nl-1).*accel_s(nt,1:nl-1,2).*rx_vector; 
        temp2 = -dm(1:nl-1).*accel_s(nt,1:nl-1,1).*ry_vector; 
        torque_inert1(nt,nl) = sum(temp1)-sum(temp2);
    clear temp1
    clear temp2
    clear rx_vector
    clear ry_vector
    end
end

Torque_totl1 = torque_inert1;
Torque_totl1(:,1)=0;
for nt = 1:num    
        for nl = 1:nl_totl
        rx_vector = X2(nt,nl:nl_totl) - x2(nt,nl);
        ry_vector = Y2(nt,nl:nl_totl) - y2(nt,nl);
        temp1 = -dm(nl:nl_totl).*accel_s(nt,nl:nl_totl,2).*rx_vector;
        temp2 = -dm(nl:nl_totl).*accel_s(nt,nl:nl_totl,1).*ry_vector;
        torque_inert2(nt,nl) = sum(temp1)-sum(temp2);
        clear temp1
        clear temp2
        clear rx_vector
        clear ry_vector
    end
end
Torque_totl2 =-torque_inert2;
Torque_totl2(:,nl_totl+1)=0;


Torque_totl = (1-s).*Torque_totl1+s.*Torque_totl2;

for nt = 1:num
    t = nt*dt;
    dcurv1(nt,:) = -a.*omega.*cos(k.*s-omega.*t);
    power(nt,:) = Torque_totl(nt,:).*dcurv1(nt,:);
    sum_power(nt) = sum(power(nt,:).*ds);
end

for nt = 1:num
    power_f(nt,:) = -dm.*accel_s(nt,:,2).*vel_int_s(nt,:,2)-dm.*accel_s(nt,:,1).*vel_int_s(nt,:,1);
    sum_power_f(nt) = sum(-power_f(nt,:));
end
figure
plot(sum_power_f,'r')
hold on
plot(sum_power,'b')
title('quxian_test')
num_sp=length(s)-2;
num_ell=40;
for n=1:num_sp
    sy(n)=0.5*(my1(n)+my1(n+1));
    sz(n)=0.5*(mz1(n)+mz1(n+1));
end

for m = 1:num_ell
    load('tri0.mat');

for n = 1:num_sp %num_sp   % straight end intersection
    
    la = abs(sy(n));
    lb = abs(sz(n));
        by(n,m) = la*cos((m-1)/num_ell*2*pi);
        bz(n,m) = lb*sin((m-1)/num_ell*2*pi);
        point((n-1)*num_ell+m, 1:3) = [s(n+1), by(n,m), bz(n,m)];
    end
end
np = length(point);
pnt_front(1,1:3) = [s(1) 0 0];
pnt_front(2:np+1,1:3) = point;
pnt_front(np+2,1:3) = [s(end) 0 0];

figure
trimesh(tri(:,:),pnt_front(:,1), pnt_front(:,2),pnt_front(:,3));
xlabel('X','FontName','Times','FontSize',20);
ylabel('Y','FontName','Times','FontSize',20);
zlabel('Z','FontName','Times','FontSize',20);
axis equal
axis tight
for n=2:num_sp
    area2(n)=0.5*(pi*sy(n-1)+pi*sz(n-1)+pi*sy(n)+pi*sz(n))*dx;
end
area2(1)=pi*(sy(1)+sz(1))/2*sqrt(dx^2+(sy(1)+sz(1))^2/4);
area2(100)=pi*(sy(99)+sz(99))/2*sqrt(dx^2+(sy(99)+sz(99))^2/4);
for nt = 1:num
    point2(:,3)= point(:,3);
    for i = 1:num_sp
        alpha_tg(nt,i) = atan((Y2(nt,i+1)-Y2(nt,i))/(X2(nt,i+1)-X2(nt,i)));
        for j=1:num_ell
            xtemp=-point((i-1)*num_ell+j,2)*sin(alpha_tg(nt,i))+x2(nt,i+1) ;
            ytemp= point((i-1)*num_ell+j,2)*cos(alpha_tg(nt,i))+y2(nt,i+1);
            point2(num_ell*(i-1)+j,1)=xtemp;
            point2(num_ell*(i-1)+j,2)=ytemp;
        end
    end

    
    pnt_front2(nt,1,1:3)=[X2(nt,1),Y2(nt,1),0];   %the accuracy of pnt_front2 is 10^-16
    pnt_front2(nt,2:np+1,1:3)=point2;
    pnt_front2(nt,np+2,1:3)=[X2(nt,100),Y2(nt,100),0];    
end
% set(0,'defaultfigurecolor','w')
% skip = 20;
% figure;
% for nt = 1:skip: 2000
%     nt
%     trimesh(tri,pnt_front2(nt,:,1), pnt_front2(nt,:,2),pnt_front2(nt,:,3));
% %         trimesh(tri,pnt(nt,:,1), pnt(nt,:,2),pnt2(nt,:,3));
% 
%     %trimesh(tri,ll(nt,:), pnt_f(nt,:,2),pnt_f(nt,:,3), value2);
%     xlabel('X','FontName','Times','FontSize',20);
%     ylabel('Y','FontName','Times','FontSize',20);
%     zlabel('Z','FontName','Times','FontSize',20);
%     hold on
%     axis equal
%     % xlim([0.7 0.9]);
%     %     xlim([-0.5 0.5]);
%     % xlim([-0.02 0.98]);
%     %     ylim([-0.2 0.2]);
%     %     zlim([-0.15 0.15])
% %         view([0 0 1])
%     %     view([1 -1 -1])
%     %     func_figure_saver(nt);
% %     hold off
% % %     pause
% end
 func_inputfiles_creation(pnt_front2);

function f = myfun(A,B,C3,nw,dt)
f = 0;
for n = 1:2*nw
    f = f + A(:,n)*B(n);
end
f = f + C3;
f = sum(abs(f))*dt;
end
