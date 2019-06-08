
for nt = 1:num
    
    %re-capture the center linew
    for nl = 2:nl_body
        nslt = (nl-2)*nbody_peri+2:(nl-1)*nbody_peri+1;
        hh(nt,nl,1:3) = mean(pnt_body_t(nt,nslt,1:3));
       
    end
     hh(nt,1,1:3)=pnt_body_t(nt,1,1:3);
     hh(nt,nl_body+1,1:3)=pnt_body_t(nt,npoint_body,1:3);
    for nl = 1:nl_body
       hhs(nt,nl,1:3) = (hh(nt,nl+1,1:3) + hh(nt,nl,1:3))/2;
    end
end


