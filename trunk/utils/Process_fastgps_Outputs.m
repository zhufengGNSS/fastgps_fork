%
% Process_fastgps_Outputs.m
%
% This program will load the output files of the fastgps software receiver
% and provide an interactive plotting interface using Octave or Matlab
%
% Copyright (c) 2008, Scott Gleason
% All rights reserved.
%
% license: GPL. For details see CD of "GNSS Applications: A Hands on Approach" in /Chapter5/fastgps/license.txt
%  

close all
clear all

disp("Loading Data Files ...")
acq_data = load('fastgps_acquisition_log_plot.dat');	
tracking_data = load('fastgps_tracking_log.dat');	
nav_data = load('fastgps_navigation_log.dat');	

%
% Interactively plot data in a loop, until Exit selected
%

% List values to plot
current_figure = 1;
channels = 99;
seconds_index = 1;
max_svs = 32;
max_channels = 12;
ch_entries = 5;     % entries per channel in the tracking data

exit_flag = 0;
while exit_flag == 0

    %
    % plotting Loop
    %

    % prompt user for packet to plot
    aaa = (input('Select: 0 - Exit, 1 - Acquisition, 2 - Tracking, 3 - Navigation -> ','s'));
    tempin = str2num(aaa);
    if tempin == 0
        break;
    else
        plot_flag1 = tempin;
    endif

    % plot
    if plot_flag1 == 1

    	% Acquisition

        figure(current_figure)
        current_figure = current_figure + 1;   % var++ and var += 1 only work in Octave

	detection_ratio = 20;	% only values of detected sats are output, i.e. depends on fastgps threshold (in parameters.h)
        data1_mag = acq_data(:,7);
        data1_ratio = acq_data(:,8);     
        data1_detected = (data1_ratio >= detection_ratio).*data1_ratio;
        data1_prn = acq_data(:,4);  

        bar(data1_prn,data1_ratio)
        bar(data1_prn,data1_detected,'g')

        xlabel('PRN')
        ylabel('Acquisition Signal Magnitude')
        title('Acquisition Results')

    elseif plot_flag1 == 2

    	% Tracking

        aaa = (input('Input Channel to Plot [0 to 11]-> ','s'));
        channel = str2num(aaa);
            
        aaa = (input('Select: 0 - Main Menu, 1 - Signal Magnitude, 2 - IQ, 3 - Doppler, 4 - Data Bits, 5 - PRN -> ','s'));
        plot_flag_tracking = str2num(aaa);

        if plot_flag_tracking == 0
            break;
        endif

        % plot selection in new figure
        figure(current_figure)
        current_figure = current_figure + 1;

        sec = tracking_data(:,seconds_index);
        if plot_flag_tracking == 1
            I = tracking_data(:,channel*ch_entries + 4);
            Q = tracking_data(:,channel*ch_entries + 5);
			signal_mag = sqrt(I.^2 + Q.^2);
            plot(sec,signal_mag,'k.')       % signal magnitude
        elseif plot_flag_tracking == 2
            I = tracking_data(:,channel*ch_entries + 4);
            Q = tracking_data(:,channel*ch_entries + 5);
            plot(I,Q,'k.')       % IQ
        elseif plot_flag_tracking == 3
            freq = tracking_data(:,channel*ch_entries + 6);
            plot(sec,freq)  % freq
        elseif plot_flag_tracking == 4
            I = tracking_data(:,channel*ch_entries + 4);
            bits = tracking_data(:,channel*ch_entries + 7);
            %plot(sec,bits)    % nav data bits, 
            plot(sec,I)    % I prompt channel 
			%axis([-inf inf -2 2])	
        elseif plot_flag_tracking == 5
            prn = tracking_data(:,channel*ch_entries + 3);
            plot(sec,prn)    % prn 
        else
            break;
        endif

    elseif plot_flag1 == 3

    	% Navigation

        aaa = (input('Select: 0 - Main Menu, 1 - Pos XYZ, 2 - Rx clock bias, 3 - Pos Lat/Lon, 4 - Height, 5 - East/North, 6 - Up, 7 - Vel XYZ, 8 - Rx clock drift, 9 - GDOP, 10 - Pseudoranges -> ','s'));
        plot_flag_nav = str2num(aaa);

        if plot_flag_nav == 0
                break;
        endif

        if plot_flag_nav == 10
            aaa = (input('Select Channel [0-11] ->','s'));
            temp_index = str2num(aaa);
	endif

        % plot selection
        sec = nav_data(:,seconds_index);
        figure(current_figure)
        current_figure = current_figure + 1;
        if plot_flag_nav == 1
            PX = nav_data(:,3);
            PY = nav_data(:,4);
            PZ = nav_data(:,5);
	    hold on
            plot(sec,PX,'k')
            plot(sec,PY,'b')
            plot(sec,PZ,'g')
        elseif plot_flag_nav == 2
            clk_bias = nav_data(:,6);
            plot(sec,clk_bias,'k')
        elseif plot_flag_nav == 3
            lat = nav_data(:,7);
            lon = nav_data(:,8);
            load mapworld.mat
            plot(world_lon,world_lat)
            plot(lon,lat,'kx')  
        elseif plot_flag_nav == 4
            hgt = nav_data(:,9);
            plot(sec,hgt,'k')
        elseif plot_flag_nav == 5
            E = nav_data(:,11);
            N = nav_data(:,12);
            plot(E,N,'kx')
        elseif plot_flag_nav == 6
            U = nav_data(:,13);
            plot(sec,U,'kx')
        elseif plot_flag_nav == 7
            VX = nav_data(:,14);
            VY = nav_data(:,15);
            VZ = nav_data(:,16);
	    hold on
            plot(sec,VX,'k')
            plot(sec,VY,'b')
            plot(sec,VZ,'g')
        elseif plot_flag_nav == 8
            clk_drift = nav_data(:,17);
            plot(sec,clk_drift,'k')
        elseif plot_flag_nav == 9
            gdop = nav_data(:,10);
            plot(sec,gdop,'k')
        elseif plot_flag_nav == 10
            PR = nav_data(:,temp_index*10 + 19);
            plot(sec,PR,'k')
        else
            break;
        endif

    else

        disp("Invalid entry, exiting")
        break

    endif % end of packet loop

endwhile  % end of while loop

% end of script


