%This program is to post-processing the eel swimming study.
clear all
% close all
clc

sub_data_input;

num = (num_end-num_start)/n_interval+1;


    %=====================================================================
    if(read_marker ~= 0 && read_marker ~= 1)
        disp('ERROR!! Wrong input for varible read_marker!!!')
      
    
    %=====================================================================
    
    elseif(read_marker == 1)
        data_body = zeros(num,npoint_body,12);
        sub_reading_markerfiles;
        
    %=====================================================================
    
    elseif(read_marker == 0)
        
        % load eel_marker,mat
        load eel_original.mat
        data_body = P0;
        [num,junk1,junk2] = size(data_body);
        
        pnt_body = data_body(:,:,1:3);
        vel_body = data_body(:,:,4:6);
        frc_body = data_body(:,:,10:12);
        
%         clear data_body data_tail p1 p2
        U_0 = 0.3;
        % obtain the velocity of fish in the lab coordinates firstly: Ux_osc, Uy_osc
        sub_velocity
        %adjust the swimming position to free swimming without incoming flow
        for nt = 1:num
            time = nt*dt*n_interval;
            pnt_body(nt,:,1) = pnt_body(nt,:,1) - U_0*time;
        end
        ux = U_0 - Ux_osc;
        uy = Uy_osc;
        U = sqrt(ux^2 + uy^2);
        % rotate & shift body 
        sub_body_rotation_shift;
%         pnt_body = pnt_body_t;
%         sub_velocity
        % obtain the centerline of mackerel with rotated positions and forces    
        sub_centerline;
        
        % obtain the position and force on each segment
        sub_force_distribution;
%         clear nslt
%         clear sts_temp
%         clear sts_int_fine
%         clear sts_int_fine_f
        
         % calculate the centerline velocity distribution
        sub_centerline_velocity;
         
         % calculate the reactivel force
        sub_reactive_force;
          
         % calculate the resistive force
         sub_resistive_force;
         clear fy_EBT
         
         %calculate the torque
         sub_torque;
         
         %save torque Torque_totl
%          clear Torque_totl1
%          clear Torque_totl2
         
         %calculate the power using torque

          sub_power_by_torque;
%          
         %calculate the power uding force
          sub_power_by_force;
             
         % calculate the details
         sub_amplitude_coefficient;   
    
         % calculate the details
         %sub_power_analysis;  
         close all
        
    end
        
