% clc
fx = sts_int_sl(:,:,1);fy = sts_int_sl(:,:,2);fz = sts_int_sl(:,:,3);
fx1 = force_inert(:,:,1);fy1 = force_inert(:,:,2);
for nt = 1:num
    ff(nt) = mean(sqrt(fx(nt,:).^2+fy(nt,:).^2));
    ffx(nt) = sum(fx(nt,:));
    ffy(nt) = sum(fy(nt,:));
    ffx1(nt) = sum(fx1(nt,:));
    ffy1(nt) = sum(fy1(nt,:));
    ffz(nt) = mean(abs(fz(nt,:)));
end
% mean(ff)
mean(ffz)
mean(ffz)/mean(ff)

% sum(ffx)
% sum(ffy)
% mean(sqrt(ffx.^2+ffy.^2))
% mean(abs(ffz))
figure
plot(ffx,'r');hold on;plot(ffx1,'k');
figure
plot(ffy,'r');hold on;plot(ffy1,'k');
