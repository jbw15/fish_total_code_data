clear all

% data = load('body/s.0000001');
data_body = load('0/s.0000001');
data_tail = load('1/s.0000001');
%tri = load('connect.dat');
load('trif.mat');
load('tri_tail.mat');


num_v = length(data_body);
num_e = length(trif);
num_v1 = length(data_tail);
num_e1 = length(tri_tail);
fid = fopen('unstruc_surface_in.dat','w');

fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%d %d\n', num_v, num_e);
fprintf(fid,'\n');
for n= 1:num_v
fprintf(fid,'%d  %f  %f  %f \n', n, data_body(n,1), data_body(n,2),data_body(n,3));
end
fprintf(fid,'\n');

for n=1:num_e
  fprintf(fid,'%d  %d  %d  %d \n', n, trif(n,1), trif(n,2),trif(n,3));  
end
fprintf(fid,'\n');
fprintf(fid, '%f %f %f \n',1, 1,1);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%d %d\n', num_v1, num_e1);
fprintf(fid,'\n');
for n= 1:num_v1
fprintf(fid,'%d  %f  %f  %f \n', n, data_tail(n,1), data_tail(n,2),data_tail(n,3));
end
fprintf(fid,'\n');

for n=1:num_e1
  fprintf(fid,'%d  %d  %d  %d \n', n, tri_tail(n,1), tri_tail(n,2),tri_tail(n,3));  
end
fprintf(fid,'\n');
fprintf(fid, '%f %f %f \n',1, 1,1);
fclose(fid);

