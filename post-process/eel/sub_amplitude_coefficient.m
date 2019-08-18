% close all
res_amp = sqrt(2*var(fy_res));
rea_amp = sqrt(2*var(fy_rea));
cfd_amp = sqrt(2*var(fy_cfd));
aa = -accel_s(:,:,2); vv = -vel_int_s(:,:,2);yh = hhs(:,:,2);
acc_amp = sqrt(2*var(aa));
vel_amp = sqrt(2*var(vv));yh_amp = sqrt(2*var(yh));
ss = 0:0.01:1;
pot = 0.5*(ss(1:end-1)+ss(2:end));
tt2 = 1:num;
for j = 2:nl_body-1
    rr3(j)=fminsearch(@(c)(sum((cfd_amp(j)*sin(2*pi/200*tt2+c)-1/3*(fy_cfd(:,j-1)'+fy_cfd(:,j)'+fy_cfd(:,j+1)')).^2)),0.01);
    rr4(j)=fminsearch(@(c)(sum((vel_amp(j)*sin(2*pi/200*tt2+c)-1/3*(vv(:,j-1)'+vv(:,j)'+vv(:,j+1)')).^2)),0.01);
    rr5(j)=fminsearch(@(c)(sum((acc_amp(j)*sin(2*pi/200*tt2+c)-1/3*(aa(:,j-1)'+aa(:,j)'+aa(:,j+1)')).^2)),0.01);
end
rr3(1  )=fminsearch(@(c)(sum((cfd_amp(1  )*sin(2*pi/200*tt2+c)-fy_cfd(:,1  )').^2)),0.01);
rr3(100)=fminsearch(@(c)(sum((cfd_amp(100)*sin(2*pi/200*tt2+c)-fy_cfd(:,100)').^2)),0.01);
rr4(1  )=fminsearch(@(c)(sum((vel_amp(1)*sin(2*pi/200*tt2+c)-vv(:,1)').^2)),0.01);
rr4(100)=fminsearch(@(c)(sum((vel_amp(100)*sin(2*pi/200*tt2+c)-vv(:,100)').^2)),0.01);
rr5(1  )=fminsearch(@(c)(sum((acc_amp(1  )*sin(2*pi/200*tt2+c)-aa(:,1  )').^2)),0.01);
rr5(100)=fminsearch(@(c)(sum((acc_amp(100)*sin(2*pi/200*tt2+c)-aa(:,100)').^2)),0.01);
rr32=rr3;rr32(1:22)=rr3(1:22)+2*pi;rr32(83:100)=rr3(83:100)-2*pi;
rr42 = rr4;rr42(1:13)=rr4(1:13)+2*pi;rr42(73:end) = rr4(73:end)-2*pi;
rr52 = rr5; rr52(1:26) = rr5(1:26)+2*pi;rr52(87:end) = rr5(87:end)-2*pi;
figure
plot(pot,rr52,'linewidth',4,'linestyle',':','color','r')
hold on
plot(pot,rr42,'linewidth',4,'linestyle','-.','color','b')
hold on
plot(pot,rr32,'linewidth',4,'linestyle','-','color','k')

xlim([0 1])
xticks([0 1])
ylim([-2*pi 2*pi])
yticks([-2*pi -pi  0 pi 2*pi])
yticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
set(gca,'FontSize',36);
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20)
set(gca,'ydir','reverse')
set(gca,'Fontsize',20,'FontName','Times')
set(gcf,'unit','centimeters','position',[8 8 15 15 ]);
gg = gcf;
print(gg,'result-figure/e-v-a-1.eps','-r600','-depsc')
hold off
% res_amp = sqrt(2*var(fy_res));
% rea_amp = sqrt(2*var(fy_rea));
% cfd_amp = sqrt(2*var(fy_cfd));
% pot=(1:100)/100;
% 
% tt2 = 1:num;
% 
% for j = 2:nl_body-1
%     rr1(j)=fminsearch(@(c)(sum((rea_amp(j)*sin(2*pi/200*tt2+c)-1/3*(fy_rea(:,j-1)'+fy_rea(:,j)'+fy_rea(:,j+1)')).^2)),0.01);
%     rr2(j)=fminsearch(@(c)(sum((res_amp(j)*sin(2*pi/200*tt2+c)-1/3*(fy_res(:,j-1)'+fy_res(:,j)'+fy_res(:,j+1)')).^2)),0.01);
%     rr3(j)=fminsearch(@(c)(sum((cfd_amp(j)*sin(2*pi/200*tt2+c)-1/3*(fy_cfd(:,j-1)'+fy_cfd(:,j)'+fy_cfd(:,j+1)')).^2)),0.01);
% end
% rr1(1  )=fminsearch(@(c)(sum((rea_amp(1  )*sin(2*pi/200*tt2+c)-fy_rea(:,1  )').^2)),0.01);
% rr1(100)=fminsearch(@(c)(sum((rea_amp(100)*sin(2*pi/200*tt2+c)-fy_rea(:,100)').^2)),0.01);
% rr2(1  )=fminsearch(@(c)(sum((res_amp(1  )*sin(2*pi/200*tt2+c)-fy_res(:,1  )').^2)),0.01);
% rr2(100)=fminsearch(@(c)(sum((res_amp(100)*sin(2*pi/200*tt2+c)-fy_res(:,100)').^2)),0.01);
% rr3(1  )=fminsearch(@(c)(sum((cfd_amp(1  )*sin(2*pi/200*tt2+c)-fy_cfd(:,1  )').^2)),0.01);
% rr3(100)=fminsearch(@(c)(sum((cfd_amp(100)*sin(2*pi/200*tt2+c)-fy_cfd(:,100)').^2)),0.01);
% 
% for j = 2:nl_body-1
%     r(j,:) = fminsearch(@(c)(sum(sum((c(1)*(fy_rea(:,j-1)+fy_rea(:,j)+fy_rea(:,j+1))+c(2)*(fy_res(:,j-1)+fy_res(:,j)+fy_res(:,j+1))-sum(fy_cfd(:,j-1:j+1),2)).^2))),[0.1,0.2]);
% end
% r(1,:) = fminsearch(@(c)(sum(sum((c(1)*(fy_rea(:,1)+fy_rea(:,2))+c(2)*(fy_res(:,1)+fy_res(:,2))-sum(fy_cfd(:,1:2),2)).^2))),[0.1,0.2]);
% r(100,:) = fminsearch(@(c)(sum(sum((c(1)*(fy_rea(:,99)+fy_rea(:,100))+c(2)*(fy_res(:,99)+fy_res(:,100))-sum(fy_cfd(:,99:100),2)).^2))),[0.1,0.2]);
% % 
% for j = 1:nl_body
%     fy_rea3(:,j) = r(j,1)*fy_rea(:,j);
%     fy_res3(:,j) = r(j,2)*fy_res(:,j);
% end
% fy_cfd3=fy_res3+fy_rea3;
% % for nl = 1:nl_body
% %     if r(nl,2)<0
% %         r(nl,2) = 0;
% %     end
% % end
% figure
% plot(pot,r(:,1),'linewidth',4,'linestyle','-.','color','b')
% hold on
% plot(pot,r(:,2),'linewidth',4,'linestyle',':','color','r')
% % % h=legend('reactive force','resistive force');
% % % set(h,'Box','off')
% % hold on
% % plot([0 1],[1 1],'Color',[0.8 0.8 0.8],'linewidth',3)
% % % hold on
% % % plot([0 1],[0 0],'Color',[0.8 0.8 0.8],'linewidth',3)
% % title('weight')
% % axis([0 1 0 10])
% % xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
% % %axis tight
% % set(gca,'FontSize',20)
% % xticks([0.01 0.505 1])
% % xticklabels({0 0.5 1})
% % xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
% % 
% figure
% plot(pot,sqrt(2*var(fy_rea3)),'linewidth',4,'linestyle','-.','color','b')
% hold on
% plot(pot,sqrt(2*var(fy_res3)),'linewidth',4,'linestyle',':','color','r')
% plot(pot,sqrt(2*var(fy_cfd)),'linewidth',4,'linestyle','-','color','k')
% % h=legend('reactive force','resistive force','CFD');
% % set(h,'Box','off')
% % xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
% % ylabel('Amplitude','FontName','Times','FontSize',25);
% title('amplitude')
% axis tight
% xlim([0 1])
% xticks([0 1])
% xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
% % 
% % set(gca,'FontSize',20)
% % 
% % rr12=rr1;rr12(1:17)=rr1(1:17)+2*pi;rr12(92:99)=rr1(92:99)-2*pi;rr12(100)=2*rr12(99)-rr12(98);%you can change the range according to your need
% % 
% % rr22=rr2;rr22(1:13)=rr2(1:13)+2*pi;rr22(72:100)=rr2(72:100)-2*pi;
% % rr32=rr3;rr32(1:22)=rr3(1:22)+2*pi;rr32(84:100)=rr3(84:100)-2*pi;
% % 
% % figure
% % plot(pot,rr12,'linewidth',4,'linestyle','-.','color','b')
% % hold on
% % plot(pot,rr22,'linewidth',4,'linestyle',':','color','r')
% % hold on
% % plot(pot,rr32,'linewidth',4,'linestyle','-','color','k')
% % % h=legend('reactive force','resistive force','CFD');
% % % set(h,'Box','off')
% % title('phase')
% % xlim([0 1])
% % xticks([0 1])
% % ylim([-2*pi 2*pi])
% % yticks([-2*pi -pi  0 pi 2*pi])
% % yticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
% % set(gca,'FontSize',36);
% % xlabel('Head<-Position->Tail','FontName','Times','FontSize',25)
% % set(gca,'ydir','reverse')
% 
% 
