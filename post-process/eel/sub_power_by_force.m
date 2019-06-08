
ds = 0.01;
for nt = 1:num
    power_hist1(nt) = sum(power_tt(nt,:).*ds);
end

 for nt = 1:num
 power_totl1(nt,:) = power_tt(nt,:).*ds;
 end
 
 for nt = 1:num
     for nl = 1:nl_totl
        power_fluid2(nt,nl) = sts_int_sl(nt,nl,2)*vel_int_s(nt,nl,2) ...
                             + sts_int_sl(nt,nl,1)*vel_int_s(nt,nl,1);
        power_inert2(nt,nl)  =   force_inert(nt,nl,2)*vel_int_s(nt,nl,2) ...
                                + force_inert(nt,nl,1)*vel_int_s(nt,nl,1);
                            
        power_fs2(nt,nl) =  power_fluid2(nt,nl)/ds;
        power_is2(nt,nl) =  power_inert2(nt,nl)/ds;
    end

 end
power_totl2 = - power_fluid2 - power_inert2 ; 
power_s2 = -power_fs2 - power_is2;
figure;
ss = 0:0.01:1;
ms = 0.5*(ss(1:end-1)+ss(2:end));
plot(ss,sum(power_positive)/2,'b','LineWidth',4)
hold on
plot(ss,sum(power_tt)/2,'-.r','LineWidth',4)
hold on
plot(ss,sum(power_sum_positive)/2,':g','LineWidth',4)
hold on
plot(ms,sum(power_s2)/2,'c--','LineWidth',4)
title('work')
% legend('W^{+}','W','W^{+}_{elasticity}','W_{fluid}')
set(gca,'ycolor','k');
ylabel('work','FontName','Times','FontSize',20);
axis tight
set(gca,'FontSize',20)
set(gca,'ycolor','k');
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);


for nt = 1:num
    power_hist2(nt) = sum(power_totl2(nt,:));
end
sum(power_hist1*dt*10)/2
sum(power_hist2*dt*10)/2

% trapz(tt,power_hist1)/2
% trapz(tt,power_hist2)/2