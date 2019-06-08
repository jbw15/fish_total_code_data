load('dm.mat')

for nt = 1:num
    for nl = 2:nl_body
        nslt = (nl-2)*nbody_peri+2:(nl-1)*nbody_peri+1;
        hh(nt,nl,1:3) = mean(pnt_body(nt,nslt,1:3));
       
    end
     hh(nt,1,1:3)=pnt_body(nt,1,1:3);
     hh(nt,nl_body+1,1:3)=pnt_body(nt,npoint_body,1:3);
    for nl = 1:nl_body
       hhs(nt,nl,1:3) = (hh(nt,nl+1,1:3) + hh(nt,nl,1:3))/2;
    end
    cx(nt)=sum(dm.*hhs(nt,:,1))/sum(dm);
    cy(nt)=sum(dm.*hhs(nt,:,2))/sum(dm);
end
Ux_osc=(cx(400)-cx(1))/(dt*400*10)
Uy_osc=(cy(400)-cy(1))/(dt*400*10)

Ux_osc1=(cx(200)-cx(1))/(dt*200*10)
Uy_osc1=(cy(200)-cy(1))/(dt*200*10) 
Ux_osc2=(cx(400)-cx(201))/(dt*200*10)
Uy_osc2=(cy(400)-cy(201))/(dt*200*10)

U1=sqrt((0.3-Ux_osc1)^2+Uy_osc1^2)
U2=sqrt((0.3-Ux_osc2)^2+Uy_osc2^2)
(sqrt((0.3-Ux_osc1)^2+Uy_osc1^2)-sqrt((0.3-Ux_osc2)^2+Uy_osc2^2))/sqrt((0.3-Ux_osc1)^2+Uy_osc1^2)

