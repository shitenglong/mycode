clc;
clear;
abs_in=0; %1 respent abs
Diameter_core_in=100;
Diameter_OC_in=200;
% shell_ratio_in=[1.2 1.5 2 2.5];%1.1:0.1:3;
shell_ratio_in=1.1:0.1:3;
% snow_radius_in=[100 200 300 500];%50:50:1000;
snow_radius_in=50:50:1000;
% BC_conc_in=0:100:5000;%0:5:5000;
BC_conc_in=0:10:1000;
result=SNICAR_Eabs_stl(Diameter_core_in,Diameter_OC_in,shell_ratio_in,...
                                   snow_radius_in,BC_conc_in,abs_in);
disp('hello,world!')
                               