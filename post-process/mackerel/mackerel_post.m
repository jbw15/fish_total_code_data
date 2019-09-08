% This program is to post-processing the mackeral swimming study. 
clear all
%clc
% close all

sub_data_input;

num = (num_end - num_start)/n_interval+1;

%==============================================================================
 if(read_marker ~= 0 && read_marker ~= 1)
     disp('ERROR!! Wrong input for varible read_marker!!!')
      %break
%==============================================================================

 elseif(read_marker == 1)
    data_body = zeros(num,npoint_body,12);
    data_tail = zeros(num,npoint_tail,12);
    sub_reading_markerfiles;
    
%==============================================================================
elseif(read_marker == 0)

    load mackerel_original_7_8.mat
    data_body = p1;
    data_tail = p2;
    [num junk1 junk2] = size(data_body);

 %   load mackeral_marker.mat
    pnt_body=data_body(:,:,1:3);
    vel_body=data_body(:,:,4:6);
    frc_body=data_body(:,:,10:12);
    
    pnt_tail=data_tail(:,:,1:3);
    vel_tail=data_tail(:,:,4:6);
    frc_tail=data_tail(:,:,10:12);
    sts_tail=data_tail(:,:,10:12);
 
    clear data_body data_tail p1 p2
    U_0 = 0.3; 
 % adjust the swimming position to free swimming without incoming flow
 % obtain the velocity of fish in the lab coordinates first: Ux_osc, Uy_osc
    sub_velocity
    for nt = 1:num
        time = nt*dt*n_interval;
        pnt_body(nt,:,1) = pnt_body(nt,:,1) - U_0*time;
        pnt_tail(nt,:,1) = pnt_tail(nt,:,1) - U_0*time;
    end
    
    load('trif.mat');
    load('trit.mat');
    % obtain the centerline of mackerel (directly from the CFD output)

    ux     = U_0 - Ux_osc;      
    uy     = Uy_osc;
    U      = sqrt(ux^2 + uy^2);

% rotate & shift body and tail
   sub_body_rotation_shift;
%     pnt_body = pnt_body_t;
%     pnt_tail = pnt_tail_t;
%     sub_velocity
%     
% refining the tail mesh for the force calculation(operate on pnt_tail and frc_tail)
    sub_refine_tail;
    [junk1 nfine_tail junk2]= size(pnt_fine_tail);

% obtain the centerline of mackerel with rotated positions and forces    
    sub_centerline_second;
    
    
% obtain the position and force on each segment
    sub_force_distribution;
    
 % calculate the centerline velocity distribution
    sub_centerline_velocity;

 % calculate the reactivel force
    sub_reactive_force;

 % calculate the resistive force
    sub_resistive_force;

    clear fy_EBT

 % calcuate the torque
    sub_torque;

 % calculate the power using torque
    sub_power_by_torque;

  % calculate the power using force
   sub_power_by_force;   
    
  % calculate the details
  sub_amplitude_coefficient;   

% close all  
 end
 
    
