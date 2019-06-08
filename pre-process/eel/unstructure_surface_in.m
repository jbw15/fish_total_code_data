clear all
clc
data = load('0/s.0000001');

% tri = load('connect.dat');
load('tri0.mat');

num_v = length(data);
num_e = length(tri);

fid = fopen('unstruc_surface_in.dat','w');

fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%d %d\n', num_v, num_e);
fprintf(fid,'\n');
for n= 1:num_v
fprintf(fid,'%d  %f  %f  %f \n', n, data(n,1), data(n,2),data(n,3));
end
fprintf(fid,'\n');

for n=1:num_e
  fprintf(fid,'%d  %d  %d  %d \n', n, tri(n,1), tri(n,2),tri(n,3));  
end
fprintf(fid,'\n');
fprintf(fid, '%f %f %f \n',1, 1,1);
fclose(fid);

