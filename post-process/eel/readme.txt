This is the total code to compute the force, torque and so on of eel:

eel_post.m is the main program, and others are subroutines;

sub_data_input.m contains the basic information of computing, such as time interval, segment of fish and so on;

sub_velocity.m computes the velocity of fish in the lab coordinate system, which will be used to computed the velocity in the incoming flow coordinate system;

sub_body_rotation_shift.m is used to rotate and shift body and tail of eel;


sub_centerline.m and sub_centerline_velocity.m compute the position, velocity and force of centerline;

sub_force_distribution.m computes the position and force of each segment;

sub_reactive_force.m and sub_resistive_force calculate the reactive force and resistive force;

sub_torque.m, sub_power_by_torque.m and sub_power_by_force calculate the torque and power by two method;

sub_amplitude_coefficient.m analysis the relationship between force computed by CFD and veloctiy and acceleration.

eel_original.mat is the data from CFD;

dm.mat , dm3.mat , sz2.mat are mass of segment, added mass of segment and height of segment;

st_fz.m is used to compute st number, kwave.m can calculate the wave speed, force_z.m is for the force in z direction;