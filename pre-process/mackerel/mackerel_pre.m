set(0,'defaultfigurecolor','w')

clear
clc

npoint=1694;
fn =  ['nlist_mackerel'];


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
pnt(:,1)=P(:,4)/(max(P(:,4))-min(P(:,4)));%X
pnt(:,2)=P(:,2)/(max(P(:,4))-min(P(:,4)));%Y
pnt(:,3)=P(:,3)/(max(P(:,4))-min(P(:,4)));%Z

% load('tri.mat');
% figure
% trimesh(tri(:,:),pnt(:,1), pnt(:,2),pnt(:,3));
% xlabel('X','FontName','Times','FontSize',20);
% ylabel('Y','FontName','Times','FontSize',20);
% zlabel('Z','FontName','Times','FontSize',20);
% axis equal
% axis tight


%=====================================================================
B=sortrows(pnt);
for i=84:1455
    pnt_f(i-83,:)=B(i,:);    
end

% figure
for j=1:49
    %         plot(pnt_f((j-1)*28+1:j*28,2),pnt_f((j-1)*28+1:j*28,3),'o')
    %         hold on
    %         plot(x_i(j,:),y_i(j,:),'ro')
    cp(j,:)=mean(pnt_f((j-1)*28+1:j*28,:));
    my1(j)=max(pnt_f((j-1)*28+1:j*28,2));
    mz1(j)=max(pnt_f((j-1)*28+1:j*28,3));
    my2(j)=min(pnt_f((j-1)*28+1:j*28,2));
    mz2(j)=min(pnt_f((j-1)*28+1:j*28,3));
    %         axis image
    %         hold off
    %         pause
end

num_ell=40;
the_i=-pi:pi/20:19*pi/20;

for j=1:49
    for i=1:28
        the(i,j)=atan2(pnt_f((j-1)*28+i,3)-cp(j,3),pnt_f((j-1)*28+i,2)-cp(j,2));
    end
    x_i(j,1:num_ell)=cp(j,1);
    y_i(j,:)=interp1(the(:,j),pnt_f((j-1)*28+1:j*28,2),the_i,'spline');
    z_i(j,:)=interp1(the(:,j),pnt_f((j-1)*28+1:j*28,3),the_i,'spline');
    
    pnt_m((j-1)*num_ell+1:j*num_ell,1)= x_i(j,:);
    pnt_m((j-1)*num_ell+1:j*num_ell,2)= y_i(j,:);
    pnt_m((j-1)*num_ell+1:j*num_ell,3)= z_i(j,:);
end

% =======================================================================%
C=B(1456:end,:);
D=sortrows(C,2);

E=D(1:135,:);

% figure
% for i=1:135
%     plot(E(i,1),E(i,3),'g*')
%     hold on
%     axis image
%     str=[num2str(i)];
%     text(E(i,1)+0.001,E(i,3),str)
%     pause
% end

ix=0;iz=0;
for i=1:135
    if E(i,3)==0
        ix=ix+1;
        xt(ix)=E(i,1); 
        pp_t(ix,:)=E(i,:);
    end
end
for i=1:135
    if E(i,2)==0&&E(i,3)~=0
        ix=ix+1;
        pp_t(ix,:)=E(i,:);
    end
end


p_t(1:19,:)=pp_t(1:19,:);
p_t(20,:)=pp_t(23,:);
p_t(21,:)=pp_t(24,:);
p_t(22,:)=pp_t(29,:);
p_t(23,:)=pp_t(30,:);
p_t(24,:)=pp_t(34,:);
p_t(25,:)=pp_t(35,:);
p_t(26,:)=pp_t(38,:);
p_t(27,:)=pp_t(39,:);

i=8;
for j=1:8
    i=i+1;
    the1=atan2(p_t(j+i,  3)-p_t(j,3),p_t(j+i,  1)-p_t(j,1));
    the2=atan2(p_t(j+i+1,3)-p_t(j,3),p_t(j+i+1,1)-p_t(j,1));
    dth1=the1:-the1/7:0;
    dth2=0:the2/7:the2;
    dth(1:8)=dth1;dth(9:15)=dth2(2:end);
    dth0=[the1,0,the2];dx0=[p_t(j+i,1),p_t(j,1),p_t(j+i+1,1)];dz0=[p_t(j+i,3),p_t(j,3),p_t(j+i+1,3)];
    
    pnt_tail0((j-1)*15+1:j*15,1)= interp1(dth0,dx0,dth,'spline');
    pnt_tail0((j-1)*15+1:j*15,2)= 0;
    pnt_tail0((j-1)*15+1:j*15,3)= interp1(dth0,dz0,dth,'spline');
end

ind0=[38 37 33 31 27 26 22 20 9 21 25 28 32 36 39];
for i=1:15
    pnt_tail0(8*15+i,:)=pp_t(ind0(i),:);
end
nt=length(pnt_tail0);
num_tail=28;
i=0;
for j = 1:2   % arc end intersection
    la = 0; %abs(sy(n));
    lb=mz1(j+47);
    for m = num_tail*3/4+1:num_tail*5/4+1
        i = i+1;
        bx(j,m) = cp(j+47,1);
        by(j,m) = la*cos((m-1)/num_tail*2*pi);
        bz(j,m) = lb*sin((m-1)/num_tail*2*pi);
        pnt_tail(i, 1:3) = [bx(j,m), by(j,m), bz(j,m)];
    end
end
pnt_tail(i+1:i+nt,:)= pnt_tail0;
pnt_tail(end,1)=1;
pnt_front(1,:)=[cp(1,1)-0.004,cp(1,2),cp(1,3)];
pnt_front(2:1961,:)=pnt_m;
pnt_front(1962,:)=[cp(end,1)+0.004,cp(end,2),cp(end,3)];

% figure
% plot3(pnt_m(:,1),pnt_m(:,2),pnt_m(:,3),'o');axis image
% hold on
% plot3(pnt_tail(:,1),pnt_tail(:,2),pnt_tail(:,3),'ro');axis image
% view([0 -1 0])

% figure;plot(pnt_tail(:,1),pnt_tail(:,3),'o');axis equal
% hold on
% for i=1:length(pnt_tail)
%     plot(pnt_tail(i,1),pnt_tail(i,3),'g*')
%     axis image
%     str=[num2str(i)];
%     text(pnt_tail(i,1)+0.001,pnt_tail(i,3),str)
%     pause
% end
%====================================================================================

% trif=tri(1:3880,:);
% trif(3881:3920,3)=1962;
% trif(3881:3920,1)=(1922:1961)';
% trif(3881:3920,2)=(1923:1962)';
% trif(3920,2)=(1922)';
load('trif.mat');
NX=11;NZ=15;
jump=0;
for j=1:NX-1
    for i=1:2*(NZ-1)
        n=2*(NZ-1)*(j-1)+i;
        tri_tail(n,1)=ceil(n/2)+jump;
        tri_tail(n,2)=tri_tail(n,1)+1;
        tri_tail(n,3)=tri_tail(n,1)+NZ+1;
        if mod(n,2)==0
            tri_tail(n,1)=tri_tail(n-1,1);
            tri_tail(n,2)=tri_tail(n-1,3);
            tri_tail(n,3)=tri_tail(n,2)-1;
        end
    end
    jump=jump+1;
end

 figure
trimesh(trif(:,:),pnt_front(:,1), pnt_front(:,2),pnt_front(:,3));

hold on
trimesh(tri_tail(:,:),pnt_tail(:,1), pnt_tail(:,2),pnt_tail(:,3));

xlabel('X','FontName','Times','FontSize',20);
ylabel('Y','FontName','Times','FontSize',20);
zlabel('Z','FontName','Times','FontSize',20);
axis equal
axis tight
% view([0 -1 0])
xt=[cp(49,1),xt,pnt_tail(160:165,1)'];

for i=2:length(xt)-1  
    for j=1:10        
        if xt(i)>=pnt_tail(j*15+15,1)&&xt(i)<pnt_tail(j*15+30,1)
            c1=-(xt(i)-pnt_tail(j*15+30,1))/(pnt_tail(j*15+30,1)-pnt_tail(j*15+15,1));
            c2= (xt(i)-pnt_tail(j*15+15,1))/(pnt_tail(j*15+30,1)-pnt_tail(j*15+15,1));
            zt(i)=c1*pnt_tail(j*15+15,3)+c2*pnt_tail(j*15+30,3);
            break
        end        
    end    
end
zt(1)=pnt_tail(30,3);

rho=1;
for i=1:length(xt)-2
    area_tail(i)=(xt(i+1)-xt(i))*(zt(i+1)+zt(i));
%     add_m(i)=rho*pi/4*(zt(i+1)+zt(i))^2*(xt(i+1)-xt(i));
    if i==10
        area_tail(i)=area_tail(i)-(xt(i+1)-xt(i))*pnt_tail(i+150,3);
%         add_m(i)=rho*pi/4*((zt(i+1)-pnt_tail(i+151,3)+zt(i))/2)^2*2*(xt(i+1)-xt(i));
    end
    if i>10
        area_tail(i)=area_tail(i)-(xt(i+1)-xt(i))*(pnt_tail(i+150,3)+pnt_tail(i+149,3));
%         add_m(i)=rho*pi/4*((zt(i+1)-pnt_tail(i+150,3)+zt(i)-pnt_tail(i+149,3))/2)^2*2*(xt(i+1)-xt(i));
    end
end
area_tail(i+1)= (zt(i+1)- pnt_tail(164,3))* (xt(i+2)-xt(i+1));

for i=1:length(xt)-2    
    add_m(i)=rho*pi/4*(zt(i+1)+zt(i))^2*(xt(i+1)-xt(i));
%     b(i) = zt(i+1)+zt(i);
    if i==10        
        add_m(i)=rho*pi/4*(zt(i)+zt(i+1)-pnt_tail(160,3))^2*(xt(i+1)-xt(i));
%         b(i)=(zt(i)+zt(i+1)-pnt_tail(160,3));
    end
    if i>10        
        add_m(i)=rho*pi/4*(zt(i)-pnt_tail(149+i,3)+zt(i+1)-pnt_tail(150+i,3))^2*(xt(i+1)-xt(i));
%         b(i)=(zt(i)-pnt_tail(149+i,3)+zt(i+1)-pnt_tail(150+i,3));
    end
end
add_m(i+1)= rho*pi/4*((zt(i+1)- pnt_tail(164,3)))^2*(xt(i+2)-xt(i+1));

% b(i+1)=(zt(i+1)- pnt_tail(164,3))/2;



for i=1:length(add_m)
      k=0;
    for j=1:length(pnt_tail)      
        if pnt_tail(j,1)>=xt(i)&&pnt_tail(j,1)<xt(i+1)
            k=k+1;
            tail_pnt(i,k)=j;
        end
    end
end

for i=2:49
    dm(i)=rho*pi/4*((my1(i)-my2(i)+my1(i-1)-my2(i-1))/2*(mz1(i)-mz2(i)+mz1(i-1)-mz2(i-1))/2)*(cp(i,1)-cp(i-1,1));%mass
%     dm3(i)=rho*pi/4*((mz1(i)-mz2(i))*(mz1(i)-mz2(i)))*(cp(i,1)-cp(i-1,1));% added mass
    dm3(i)=rho*pi/4*((mz1(i)-mz2(i)+mz1(i-1)-mz2(i-1))/2)^2*(cp(i,1)-cp(i-1,1));% added mass
%     a(i) = (mz1(i)-mz2(i)+mz1(i-1)-mz2(i-1))/2;
end
%==================================================================================
dm(1)=1/3*0.004*pi/4*(my1(1)-my2(1))*(mz1(1)-mz2(1))*rho;         %dm(2)*0.2;
dm3(1)=pi*((mz1(1)-mz2(1))/2)^2*0.004*rho; %dm3(2)*0.2;
% a(1) = mz1(1)-mz2(1);
j=length(add_m);
dm(i+1:i+j)=0;                  %dm(1);
% dm(i+1:i+j)=dm(1);
dm3(i+1:i+j)=add_m;
% a(i+1:i+j)=b;

% for i=50:58
%     if i==50
%         dm3(i)=rho*pi* pnt_tail0((i-49)*15,3)^2*(xt(1)-cp(end,1));
%     else
%         dm3(i)=rho*pi* pnt_tail0((i-49)*15,3)^2*(xt(i-50+1)-xt(i-50));
%     end
%     sz2(i)=pnt_tail0((i-49)*15,3)*2;
% end
% sz2(1:49)=mz1-mz2;
clear B C D E
%====================================================
m=sum(dm);

ss=[cp(1,1)-0.004,cp(1:49,1)',xt(2:end)];
dx=[0,diff(ss)];
for i=2:49
    area_front(i)=(ss(i+1)-ss(i))*(mz1(i)-mz2(i)+mz1(i-1)-mz2(i-1))/2;
end
area_front(1)=0.4*area_front(2);
area=[area_front,area_tail];

clear c1 c2
lamd=1.0;                                                                                                                                      
a0=0.02*50;a1=-0.08*40;a2=0.16*35;% case: mackerel
k=2*pi/lamd;w=2*pi;
dt=5e-4;
Numt=round(2*pi/w/dt);
Numx=length(ss);

dt = 5e-4; t = dt:dt:1;num = length(t);
ds = dx; s = ss; nl_totl = length(s) - 1;
omega = 2*pi; lambda =1; k = 2*pi/lambda;
a0 = 1; a1 = -3.2; a2 = 5.6; a = a0 + a1*s + a2*s.^2;
nw=3;
for nt = 1:num
    t = nt*dt;
    curv(nt,:) = a.*sin(k*s-omega*t);
    for nl = 2:nl_totl+1
        xx = [s(1:nl)];the = [curv(nt,1:nl)];
        theta(nt,nl) = trapz(xx,the);
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
    cx1(nt)=sum(dm.*X1(nt,:))/m; % The position of COM is at origin of coordinates after translation
    cy1(nt)=sum(dm.*Y1(nt,:))/m; % momentum is conserved after shifting fish.
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
%     A(nt,1:nw) = C1(nt)*cosw(nt,:);
%     A(nt,nw+1:2*nw) = C1(nt)*sinw(nt,:);

end 
b0 = [0.3293   -0.0000  0  0.9027   -0.0000 0] ;
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
   % x1(nt,:)=x1(nt,:)-cx1(nt);     
   % y1(nt,:)=y1(nt,:)-cy1(nt);    
    x2(nt,:) = x1(nt,:)*cos(rot(nt)) - y1(nt,:)*sin(rot(nt)); %The position of the center of mass does not change when fish rotates around the center of mass
    y2(nt,:) = y1(nt,:)*cos(rot(nt)) + x1(nt,:)*sin(rot(nt));   
    X2(nt,:) = X1(nt,:)*cos(rot(nt)) - Y1(nt,:)*sin(rot(nt));
    Y2(nt,:) = Y1(nt,:)*cos(rot(nt)) + X1(nt,:)*sin(rot(nt));    
    cx2(nt)=sum(dm.*X2(nt,:))/m;
    cy2(nt)=sum(dm.*Y2(nt,:))/m;
end

for nt = 1:num
    X2(nt,:) = X2(nt,:) - cx2(nt); 
    Y2(nt,:) = Y2(nt,:) - cy2(nt);
    x2(nt,:) = x2(nt,:) - cx2(nt);
    y2(nt,:) = y2(nt,:) - cy2(nt);
    %Y2(nt,:) = 0.5*(y2(nt,1:nl_totl)+y2(nt,2:nl_totl+1));
    ic2(nt,:)=dm.*(X2(nt,:).^2+Y2(nt,:).^2); % The position of COM is at the origin 
    Icom2(nt)=sum(ic2(nt,:));
end

for i=2:Numt-1
    DI(i)=(Icom2(i+1)-Icom2(i-1))/(2*dt);
end
DI(1)=2*DI(2)-DI(3);
DI(Numt)=2*DI(Numt-1)-DI(Numt-2);


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
nl_connect = nl_totl/2;

% Torque_totl(:,1:nl_connect) = Torque_totl1(:,1:nl_connect);
% Torque_totl(:,nl_connect+1:nl_totl+1) = -Torque_totl2(:,nl_connect+1:nl_totl+1);

Torque_totl = (1-s).*Torque_totl1+s.*Torque_totl2;
%Torque_totl = Torque_totl.*ds;
% tt = linspace(0,1,num);
% xx = linspace(0,1,nl_totl+1);
% [xxx ttt] = meshgrid(xx,tt);
% 
% torque = -0.08359.*((0.0208333 + (1.28245e-16).*xxx - 0.0208333.*cos(12.5664.*xxx)).*sin(6.28319.*ttt) + cos(6.28319.*ttt).*((-5.55112e-17) + 0.261799.*xxx -...
%         0.785398.*xxx.^2 +0.523599.*xxx.^3 - 0.0208333.*sin(12.5664.*xxx)));


for nt = 1:num
    t = nt*dt;
    dcurv1(nt,:) = -a.*omega.*cos(k.*s-omega.*t);
%     if t<0.25
%         dcurv1(nt,:) = -omega_2.*a.*cos(k*s-omega_2*t);
%     elseif (t>=0.25)&&(t<=0.75)
%         dcurv1(nt,:) = -omega_1.*a.*cos(k*s-omega_1*t-pi);
%     else
%         dcurv1(nt,:) = omega_2.*a.*cos(k*s-omega_2*t);
%     end
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

np=length(pnt_front);
for nt = 1:num
    point2(:,3)= pnt_front(2:np-1,3);
    for i = 1:49
        alpha_tg(nt,i) = atan((Y2(nt,i+1)-Y2(nt,i))/(X2(nt,i+1)-X2(nt,i)));
        for j=1:num_ell
            xtemp=-pnt_front((i-1)*num_ell+j+1,2)*sin(alpha_tg(nt,i))+x2(nt,i+1) ;
            ytemp= pnt_front((i-1)*num_ell+j+1,2)*cos(alpha_tg(nt,i))+y2(nt,i+1);
            point2(num_ell*(i-1)+j,1)=xtemp;
            point2(num_ell*(i-1)+j,2)=ytemp;
        end
    end
    pnt_front2(nt,1,1:3)=[X2(nt,1),Y2(nt,1),cp(1,3)];
    pnt_front2(nt,2:np-1,1:3)=point2;
    pnt_front2(nt,np,1:3)=[X2(nt,50),Y2(nt,50),cp(49,3)];% accuracy is 10^-16
end
%===========================Tail========================================
for nt = 1:num
    pnt_tail2(nt,:,3)= pnt_tail(:,3);
    pnt_tail2(nt,1:15,1)=mean(pnt_front2(nt,47*40+2:48*40+1,1));
    pnt_tail2(nt,1:15,2)=mean(pnt_front2(nt,47*40+2:48*40+1,2));
    pnt_tail2(nt,16:30,1)=mean(pnt_front2(nt,48*40+2:49*40+1,1));
    pnt_tail2(nt,16:30,2)=mean(pnt_front2(nt,48*40+2:49*40+1,2));
    pnt_tail2(nt,31:end,1)=interp1(ss,x2(nt,:),pnt_tail(31:end,1),'spline');
    pnt_tail2(nt,31:end,2)=interp1(ss,y2(nt,:),pnt_tail(31:end,1),'spline');
end



func_inputfiles_creation0(pnt_front2);
func_inputfiles_creation1(pnt_tail2);

function f = myfun(A,B,C3,nw,dt)
f = 0;
for n = 1:2*nw
    f = f + A(:,n)*B(n);
end
f = f + C3;
f = sum(abs(f))*dt;
end
