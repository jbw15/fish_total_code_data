
    for nt = 1:num
        
        % re-capture the center line
        for nl=2:nl_body
            nslt = (nl-2)*nbody_peri+2:(nl-1)*nbody_peri+1;
            hh(nt,nl,1:3) = mean(pnt_body_t(nt,nslt, 1:3));
        end
            hh(nt,1,1:3) = pnt_body_t(nt,1, 1:3);

        % tail before forking
        for nl=nl_body+1:nl_body+nl_tail_1-3
            nslt = (nl-nl_body+1)*ntail_vtcl+8;
            hh(nt,nl,1:3) = pnt_tail_t(nt,nslt, 1:3);
        end
            nl = nl_body+nl_tail_1-2;
            nslt = (nl-nl_body+1)*ntail_vtcl+9;
            hh(nt,nl,1:3) = pnt_tail_t(nt,nslt, 1:3);
        % tail after forking
         for nl = nl_body+nl_tail_1-1:nl_totl
            hh(nt,nl,1:3) =  pnt_tail_t(nt,nl+100, 1:3);
         end
           
      for nl = 1:nl_totl-1
         hhs(nt,nl,1:3) = (hh(nt,nl+1,1:3) + hh(nt,nl,1:3))/2;
      end
    end

