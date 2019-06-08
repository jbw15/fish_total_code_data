a1=max(hhs(1:200,end,2));
a2=min(hhs(1:200,end,2));
st1 = (a1-a2)/U1;
a1=max(hhs(201:400,end,2));
a2=min(hhs(201:400,end,2));
st2 = (a1-a2)/U2;
st = (st1+st2)/2
% fz = sts_int_sl(:,:,3);
% for nt = 1:num
% ffz(nt) = sum(fz(nt,:));
% end
% mean(abs(ffz))