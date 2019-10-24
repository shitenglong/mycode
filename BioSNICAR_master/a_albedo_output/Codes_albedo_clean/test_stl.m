clear;clc;
% FILE NAMES CONTAINING MIE PARAMETERS FOR ALL AEROSOL SPECIES:
    % (ideally, these files should exist in all 'band' directories)
    fl_sot1  = 'mie_sot_ChC90_dns_1317.nc';
    fl_sot2  = 'miecot_slfsot_ChC90_dns_1317.nc';
    fl_dst1  = 'aer_dst_bln_20060904_01.nc';
    fl_dst2  = 'aer_dst_bln_20060904_02.nc';
    fl_dst3  = 'aer_dst_bln_20060904_03.nc';
    fl_dst4  = 'aer_dst_bln_20060904_04.nc';
    fl_ash1  = 'volc_ash_mtsthelens_20081011.nc';
    fl_bio1  = 'biological_1.nc'; % Biological impurity 1 (30um diameter, 1.5%chll a, 10% each 1 % 2 carotenoids)
    fl_bio2  = 'biological_2.nc'; % Biological impurity 2 (30um diameter, 1.5%chll a, 5% each 1 % 2 carotenoids)
    fl_bio3  = 'biological_3.nc'; % Biological impurity 3 (30um diameter, 1.5%chll a, 1% each 1 % 2 carotenoids)
    fl_bio4  = 'biological_4.nc'; % Biological impurity 4 (30um diameter, 1.5%chll a only)
    fl_bio5  = 'biological_5.nc'; % Biological impurity 5 (10um diameter, pigs as per bio2)
    fl_bio6  = 'biological_6.nc'; % Biological impurity 6 (50um diameter, pigs as per bio2)
    fl_bio7  = 'biological_7.nc'; % Biological impurity 6 (20um diameter, pigs as per bio2)
    fl_water1 = 'water_segelstein_20.nc'; %water type 1 (0.2mm spheres)
    
    
    % COSINE OF SOLAR ZENITH ANGLE FOR DIRECT-BEAM RT
%     coszen=0.5; %需要修改
  %  coszen=coszen_in; %需要修改 @@@@@@@@@@@@@@@@@@@@@天顶角的余弦值
  
    % REFLECTANCE OF SURFACE UNDERLYING SNOW:
    % (value applied to all wavelengths.  user can also specify
    % spectrally-dependent ground albedo below)
%     R_sfc_all_wvl = 0.25; %需要修改
   R_sfc_all_wvl = 0.35; %需要修改
  
    % RADIATIVE TRANSFER CONFIGURATION:
    BND_TYP  = 1;        % 1= 470 spectral bands, 2= 5 bands, 3= 3 bands
    DIRECT   =1;% DIRECT_IN;        % 1= Direct-beam incident flux, 0= Diffuse incident flux
    APRX_TYP = 3;        % 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    DELTA    = 1;        % 1= Apply Delta approximation, 0= No delta
    
    
   

wrkdir_aero='E:\2018ACP_pw\BioSNICAR_master\Impurity_RIs\';
fl_in1 = strcat(wrkdir_aero,fl_sot1);
fl_in2 = strcat(wrkdir_aero,fl_sot2);
fl_in3 = strcat(wrkdir_aero,fl_dst1);
fl_in4 = strcat(wrkdir_aero,fl_dst2);
fl_in5 = strcat(wrkdir_aero,fl_dst3);
fl_in6 = strcat(wrkdir_aero,fl_dst4);
fl_in7 = strcat(wrkdir_aero,fl_ash1);
fl_in8 = strcat(wrkdir_aero,fl_bio1);
fl_in9 = strcat(wrkdir_aero,fl_bio2);
fl_in10 = strcat(wrkdir_aero,fl_bio3);
fl_in11 = strcat(wrkdir_aero,fl_bio4);
fl_in12 = strcat(wrkdir_aero,fl_bio5);
fl_in13 = strcat(wrkdir_aero,fl_bio6);
fl_in14 = strcat(wrkdir_aero,fl_water1);
fl_in15 = strcat(wrkdir_aero,fl_bio7);

omega_aer(:,1)       = ncread(fl_in1,'ss_alb');
g_aer(:,1)           = ncread(fl_in1,'asm_prm');
ext_cff_mss_aer(:,1) = ncread(fl_in1,'ext_cff_mss');
  
omega_aer(:,2)       = ncread(fl_in2,'ss_alb');
g_aer(:,2)           = ncread(fl_in2,'asm_prm');
ext_cff_mss_aer(:,2) = ncread(fl_in2,'ext_cff_mss_cor');
%ext_cff_mss_aer(:,2)= ncread(fl_in2,'ext_cff_mss');

omega_aer(:,3)       = ncread(fl_in3,'ss_alb');
g_aer(:,3)           = ncread(fl_in3,'asm_prm');
ext_cff_mss_aer(:,3) = ncread(fl_in3,'ext_cff_mss');

omega_aer(:,4)       = ncread(fl_in4,'ss_alb');
g_aer(:,4)           = ncread(fl_in4,'asm_prm');
ext_cff_mss_aer(:,4) = ncread(fl_in4,'ext_cff_mss');

omega_aer(:,5)       = ncread(fl_in5,'ss_alb');
g_aer(:,5)           = ncread(fl_in5,'asm_prm');
ext_cff_mss_aer(:,5) = ncread(fl_in5,'ext_cff_mss');

omega_aer(:,6)       = ncread(fl_in6,'ss_alb');
g_aer(:,6)           = ncread(fl_in6,'asm_prm');
ext_cff_mss_aer(:,6) = ncread(fl_in6,'ext_cff_mss');

omega_aer(:,7)       = ncread(fl_in7,'ss_alb');
g_aer(:,7)           = ncread(fl_in7,'asm_prm');
ext_cff_mss_aer(:,7) = ncread(fl_in7,'ext_cff_mss');

omega_aer(:,8)       = ncread(fl_in8,'ss_alb');
g_aer(:,8)           = ncread(fl_in8,'asm_prm');
ext_cff_mss_aer(:,8) = ncread(fl_in8,'ext_cff_mss');

omega_aer(:,9)       = ncread(fl_in9,'ss_alb');
g_aer(:,9)           = ncread(fl_in9,'asm_prm');
ext_cff_mss_aer(:,9) = ncread(fl_in9,'ext_cff_mss');

omega_aer(:,10)       = ncread(fl_in10,'ss_alb');
g_aer(:,10)           = ncread(fl_in10,'asm_prm');
ext_cff_mss_aer(:,10) = ncread(fl_in10,'ext_cff_mss');

omega_aer(:,11)       = ncread(fl_in11,'ss_alb');
g_aer(:,11)           = ncread(fl_in11,'asm_prm');
ext_cff_mss_aer(:,11) = ncread(fl_in11,'ext_cff_mss');

omega_aer(:,12)       = ncread(fl_in12,'ss_alb');
g_aer(:,12)           = ncread(fl_in12,'asm_prm');
ext_cff_mss_aer(:,12) = ncread(fl_in12,'ext_cff_mss');

omega_aer(:,13)       = ncread(fl_in13,'ss_alb');
g_aer(:,13)           = ncread(fl_in13,'asm_prm');
ext_cff_mss_aer(:,13) = ncread(fl_in13,'ext_cff_mss');

omega_aer(:,14)       = ncread(fl_in14,'ss_alb');
g_aer(:,14)           = ncread(fl_in14,'asm_prm');
ext_cff_mss_aer(:,14) = ncread(fl_in14,'ext_cff_mss');

omega_aer(:,15)       = ncread(fl_in15,'ss_alb');
g_aer(:,15)           = ncread(fl_in15,'asm_prm');
ext_cff_mss_aer(:,15) = ncread(fl_in15,'ext_cff_mss');

% Set aerosol concentration matrix:
mss_cnc_aer(1:nbr_lyr,1) = mss_cnc_sot1;
mss_cnc_aer(1:nbr_lyr,2) = mss_cnc_sot2;
mss_cnc_aer(1:nbr_lyr,3) = mss_cnc_dst1;
mss_cnc_aer(1:nbr_lyr,4) = mss_cnc_dst2;
mss_cnc_aer(1:nbr_lyr,5) = mss_cnc_dst3;
mss_cnc_aer(1:nbr_lyr,6) = mss_cnc_dst4;
mss_cnc_aer(1:nbr_lyr,7) = mss_cnc_ash1;
mss_cnc_aer(1:nbr_lyr,8) = mss_cnc_bio1;
mss_cnc_aer(1:nbr_lyr,9) = mss_cnc_bio2;
mss_cnc_aer(1:nbr_lyr,10) = mss_cnc_bio3;
mss_cnc_aer(1:nbr_lyr,11) = mss_cnc_bio4;
mss_cnc_aer(1:nbr_lyr,12) = mss_cnc_bio5;
mss_cnc_aer(1:nbr_lyr,13) = mss_cnc_bio6;
mss_cnc_aer(1:nbr_lyr,14) = mss_cnc_water1;
mss_cnc_aer(1:nbr_lyr,15) = mss_cnc_bio7;

% convert to units of kg/kg:
mss_cnc_aer = mss_cnc_aer.*10^-9;

