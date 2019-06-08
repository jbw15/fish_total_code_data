for nt=1:num
    for nl = 1:8  
        for j = 1:ntail_vtcl
            
            ii = (nl+1)*ntail_vtcl+j;
            jj = (nl+2)*ntail_vtcl+j;
            point_temp = (pnt_tail_t(nt,ii,1:3) + pnt_tail_t(nt,jj,1:3))/2.0;
            vel_temp   = (vel_tail_t(nt,ii,1:2) + vel_tail_t(nt,jj,1:2))/2.0;
            sts_temp   = (sts_tail_t(nt,ii,1:3) + sts_tail_t(nt,jj,1:3))/2.0;
            frc_temp   = (frc_tail_t(nt,ii,1:3) + frc_tail_t(nt,jj,1:3))/4.0;
            
            kk1 = (nl-1)*2*ntail_vtcl+j;
            kk2 = (nl*2-1)*ntail_vtcl+j;
            
            pnt_fine_tail(nt,kk1,1:3) =  pnt_tail_t(nt,ii,1:3);
            pnt_fine_tail(nt,kk2,1:3) =  point_temp;
            
            vel_fine_tail(nt,kk1,1:2) =  vel_tail_t(nt,ii,1:2);
            vel_fine_tail(nt,kk2,1:2) =  vel_temp;
            
            frc_fine_tail(nt,kk1,1:3) = frc_tail_t(nt,ii,1:3)/2;
            frc_fine_tail(nt,kk2,1:3) = frc_temp;
            
            sts_fine_tail(nt,kk1,1:3) =  sts_tail_t(nt,ii,1:3);
            sts_fine_tail(nt,kk2,1:3) =  sts_temp;
            
        end 
    end
    nl = 9;
    for j= 1:ntail_vtcl
        ii = (nl+1)*ntail_vtcl+j;
        kk1 = (nl-1)*2*ntail_vtcl+j;
        pnt_fine_tail(nt,kk1,1:3) =  pnt_tail_t(nt,ii,1:3);
        frc_fine_tail(nt,kk1,1:3) =  frc_tail_t(nt,jj,1:3)/2;
        sts_fine_tail(nt,kk1,1:3) =  sts_tail_t(nt,ii,1:3);
        vel_fine_tail(nt,kk1,1:2) =  vel_tail_t(nt,ii,1:2);
    end
end
% 
% figure;plot(pnt_fine_tail(nt,:,1),pnt_fine_tail(nt,:,3),'o');axis equal
% hold on
% for i=1:length(pnt_fine_tail)
%     plot(pnt_fine_tail(nt,i,1),pnt_fine_tail(nt,i,3),'g*')
%     axis image
%     str=[num2str(i)];
%     text(pnt_fine_tail(nt,i,1)+0.001,pnt_fine_tail(nt,i,3),str)
%     pause
% end

