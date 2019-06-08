for nt = 2:num-1
     if nt>1 && nt < num
        for j=1:3
            vel_int_s(nt,:,j) = (hhs(nt+1,:,j)-hhs(nt-1,:,j))/(2*dt*n_interval);
        end
     end
 end

 for j = 1:3
     vel_int_s(1,:,j) = 2*vel_int_s(2,:,j) - vel_int_s(3,:,j);
     vel_int_s(num,:,j) = 2*vel_int_s(num-1,:,j) - vel_int_s(num-2,:,j);
 end

% figure();
% imagesc(pot,xt,vel_int_s(:,:,2));
% axis xy;
% colorbar;
% xlabel('Head<-Position->Tail','FontName','Times');
% ylabel('Time','FontName','Times');
% title('lateral velocity distribution')
% set(gca,'FontSize',20);
% colormap('jet');
% xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
% ylabel('\fontsize{20}\fontname{Times new roman}Time')

