filename = 'input_file.dat';
fid = fopen(filename,'r');
disp('=====================================================')
fprintf('Reading datas from file %s ...\n', filename);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
read_marker = sscanf(line, '%d');
fprintf(1,'read_marker: %d \n',read_marker);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
temp_input = sscanf(line, '%d');
num_start = temp_input(1);
num_end = temp_input(2);
fprintf('num_start, num_end: [%d %d] \n', num_start, num_end);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
npoint_body = sscanf(line, '%d');
fprintf(1,'npoint_body: %d \n',npoint_body);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
n_interval = sscanf(line, '%d');
fprintf(1,'n_interval: %d \n',n_interval);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
nbody_peri = sscanf(line, '%d');
fprintf(1,'nbody_peri: %d \n',nbody_peri);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
nl_body = sscanf(line, '%d');
fprintf('nl_body: %d \n', nl_body);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
nl_totl = sscanf(line, '%d');
fprintf(1,'nl_totl: %d \n',nl_totl);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
temp_input = sscanf(line,'%f');
Ux_osc = temp_input(1);
Uy_osc = temp_input(2);
fprintf(1,'Ux_osc, Uy_osc: %f %f \n',Ux_osc,Uy_osc);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
U_0 = sscanf(line, '%f');
fprintf(1,'U_0: %8.4f \n',U_0);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
dt = sscanf(line, '%f');
fprintf(1,'dt: %5.2e \n',dt);
%========================================================
line = fgetl(fid);
line = fgetl(fid);
temp_input = sscanf(line, '%d');
refine1 = temp_input(1);
refine2 = temp_input(2);
fprintf('refine1, refine2: %d %d \n', refine1, refine2);
%========================================================
fclose(fid);
disp('finish reading data from the file...')
disp('========================================================')
