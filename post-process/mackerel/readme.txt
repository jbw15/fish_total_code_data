﻿This is the total code to compute the force, torque and so on of mackerel:

mackerel_post.m is the main program, and others are subroutines;

sub_data_input.m contains the basic information of compute, such as time interval, segment of fish and so on;

sub_velocity.m computes the velocity of fish in the lab coordinate system, which will be used to computed the velocity in the incoming flow coordinate system;

sub_body_rotation_shift.m is used to rotate and shift body and tail of mackerel;

sub_refine_tail.m is used to refine the tail mesh for the force calculation;

sub_centerline_second.m and sub_centerline_velocity.m compute the position, velocity and force of centerline;

sub_force_distribution.m computes the position and force of each segment;

sub_reactive_force.m and sub_resistive_force calculate the reactive force and resistive force;

sub_torque.m, sub_power_by_torque.m and sub_power_by_force calculate the torque and power by two method;

sub_amplitude_coefficient.m analysis the relationship between force computed by CFD and velocity and acceleration.

mackerel_original.mat and mackerel_original_7_8.mat are the data from CFD;(different two cycles.)

dm.mat , dm3_1.mat , ss.mat, area.mat are mass of segment, added mass of segment arc length of segment and area of segment.

st_fz.m is used to compute st number, kwave.m can calculate the wave speed, force_z.m is for the force in z direction;