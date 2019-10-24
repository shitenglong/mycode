clc;
clear;
abs_in=0; % 1respent abs   0 nonabs
Diameter_core_in=100;
Diameter_OC_in=200;
shell_ratio_in=[1.2 1.5 2 2.5];
snow_radius_in=[100 200 300 500];
% BC_conc_in=0:100:5000;
BC_conc_in=0:10:1000;
result=SNICAR_Eabs_stl_v(Diameter_core_in,Diameter_OC_in,shell_ratio_in,...
                                   snow_radius_in,BC_conc_in,abs_in);
disp('hello,world!')
                               