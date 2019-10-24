clc;clear;
%eval('cd E:\2018ACP_pi\ACP\');
coszen=0.5;
bc_oc=csvread('E:\2018ACP_pw\ACP\shitlwork\input_time\site_time.csv',1,0);
oc_bc_ratio=bc_oc(:,5)./bc_oc(:,4);
index=find(oc_bc_ratio > 1 & oc_bc_ratio < 30);
number_sample=size(index,1);
site_select=bc_oc(index,1);
oc_bc_select=oc_bc_ratio(index);
bc_select=bc_oc(index,4);
oc_select=bc_oc(index,5);
lon_select=bc_oc(index,2);
lat_select=bc_oc(index,3);
% bc_select=[0 ;bc_select];
% oc_select=[0; oc_select];
oc_bc_select=[1; oc_bc_select];%252
%  ======================================
m_BC=complex(1.85,0.71);
% m_OC=complex(1.55,0.001);
d_bc=1.8;%  g cm-3
d_oc=1.2;
wave_length=305:10:4995;
k_OC_all=load('E:\2018ACP_pw\ACP\shitlwork\Mie_model\K_OC\300_750_K_OCk_brown_0.06_AAE_brown_-6_Diameter_cores_200nm.txt');
k_OC_all=[k_OC_all(6:10:451);zeros(1,425)'];
Diameter_core=100;%[100 200 300 400 500]; %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@xiugai
% shell_ratio=[1.5 2 3 4];
size_wave=size(wave_length);
snow_radius=100;%[100 200 500 1000];
%    循环开始
for d=1:1;
   for r=1:1;
     for s=90:number_sample;
         clearvars -except oc_bc_ratio number_sample site_select oc_bc_select bc_select oc_select lon_select lat_select ...
         m_BC d_bc d_oc wave_length k_OC_all Diameter_core size_wave snow_radius d r s  
      % ===========================求太阳高度角
      if (site_select(s) >= 1 & site_select(s) <= 46)
               NF = 2010;
      elseif (site_select(s) >= 47 & site_select(s) <= 84)
               NF = 2012;
      else
               NF = 2014;
      end
% NF = 2010; %年份 
Y = 1;  %月
% R = 25;  %日
R = 1:1:31;%[17 17 19 19 25 25 25];
D = lon_select(s);  %观测处经度的度值
E =0; %0.4*60;  %观测处经度的分度值
phi =lat_select(s);%观测地的纬度
%S=16;  %观测时的时值
%F=32;  %观测时的分值
S =4:1:22;
F=0;%F = [20 30 56 11 41 47 32];
% 计算t
N0 = 79.6764+0.2422*(NF-1985)-floor((NF-1985)/4);
%积日N
A = NF/4;
B = A-floor(A);
C = 32.8;
if Y <= 2
    C = 30.6;
end
if (B == 0) && (Y > 2)
    C = 31.8;
end
        for ii=1:31; 
            % ========
 N = floor(30.6*Y-C+0.5)+R(ii);
%积日修正值dN
L = (D+E/60)/15;        %经度订正
W = (S+F)/60;           %时刻订正
dN = (W-L)/24;
t = N-N0+dN;
%t = (floor(30.6*Y-C+0.5)+R)+((S+F/60)-(D+E/60)*15)/24-(79.6764+0.2422*(NF-1985)-floor((NF-1985)/4));
%  日角
theta = 2*pi*t/365.2422;
theta = theta*(pi/180);
%  太阳赤纬
delta = 0.3723+23.2567*sin(theta)+0.1149*sin(2*theta)-0.1712*sin(3*theta)-0.758*cos(theta)+0.3656*cos(2*theta)+0.0201*cos(3*theta);
%  时差
Eq = 0.0028-1.9857*sin(theta)+9.9059*sin(2*theta)-7.0924*cos(theta)-0.6882*cos(2*theta);
% 太阳时角
LC = 4*(D-120);
TT = (S+F/60+LC/60+Eq/60);
omega = (TT-12)*15;
% 转换为弧度制
delta = delta*(pi/180);
phi = phi*(pi/180);
omega = omega*(pi/180);
%  太阳高度角
hs = asind(sin(phi)*sin(delta)+cos(phi)*cos(delta).*cos(omega));
%   完毕
indexs_hour=find(hs > 0);
hs_select=hs(indexs_hour);
S_select=S(indexs_hour);
numbers_hour=size(indexs_hour);
coszen=sin(hs_select./180*pi);
% ======================================
        shell_ratio=(oc_bc_select(s)*d_bc/d_oc+1)^(1/3);
        Diameter_OC=Diameter_core(d)*(oc_bc_select(s)*d_bc/d_oc)^(1/3);
        Diameter_shell=shell_ratio*Diameter_core(d);
        d_shell=(Diameter_core(d)/Diameter_shell)^3*d_bc+(1-(Diameter_core(d)/Diameter_shell)^3)*d_oc;
        BC_conc_factor=d_shell/d_bc*(Diameter_shell/Diameter_core(d))^3;
        
          for i=1:size_wave(2);
            Dx_wave_core=Diameter_core(d)*pi/wave_length(i);
            Dx_wave_shell=Diameter_shell*pi/wave_length(i);
            Dx_wave_OC=Diameter_OC*pi/wave_length(i);%@@@@@@@@@@@@@@@@@@@@@@@@
            m_OC=complex(1.55,k_OC_all(i));          %@@@@@@@@@@@@@
           % Cabs_OC_difference(i)=result_shell(3)*pi*poier(Diameter_shell,2)/4-result_core(3)*pi*poier(Diameter_core,2)/4;
            result_core_BC=Mie(m_BC,Dx_wave_core);
            Cabs_core_BC(i)=result_core_BC(3)*pi*power(Diameter_core(d),2)/4/(pi*Diameter_core(d)^3/6.0*d_bc);
            result_core_BC_shell_OC=Miecoated(m_BC,m_OC,Dx_wave_core,Dx_wave_shell,1);
            Cabs_core_BC_shell_OC(i)=result_core_BC_shell_OC(3)*pi*power(Diameter_shell,2)/4/(pi*Diameter_shell^3/6.0*d_shell);
            Eabs(i)=BC_conc_factor*Cabs_core_BC_shell_OC(i)/Cabs_core_BC(i);
            data_bc_core_spec_in(i,1)=result_core_BC(2)/result_core_BC(1);
            data_bc_core_spec_in(i,2)=result_core_BC(5);
            data_bc_core_spec_in(i,3)=result_core_BC(1)*pi*power(Diameter_core(d),2)/4/(pi*Diameter_core(d)^3/6.0*d_bc)*10^6;
            %result_core_BC_shell_nonabs=Miecoated(m_BC,m_nonabs,Dx_iave_core,Dx_iave_shell,1);
            %Cabs_core_BC_shell_nonabs(i)=result_core_BC_shell_nonabs(3)*pi*poier(Diameter_shell,2)/4;
            data_bc_shell_spec_in(i,1)=result_core_BC_shell_OC(2)/result_core_BC_shell_OC(1);
            data_bc_shell_spec_in(i,2)=result_core_BC_shell_OC(5);
            data_bc_shell_spec_in(i,3)=result_core_BC_shell_OC(1)*pi*power(Diameter_shell,2)/4/(pi*Diameter_shell^3/6.0*d_shell)*10^6;
            %+++++++++++++++++++++++++++++++++++++++++++++++++oc@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            result_core_OC=Mie(m_OC,Dx_wave_OC);
            data_bc_OC_spec_in(i,1)=result_core_OC(2)/result_core_OC(1);
            data_bc_OC_spec_in(i,2)=result_core_OC(5);
            data_bc_OC_spec_in(i,3)=result_core_OC(1)*pi*power(Diameter_OC,2)/4/(pi*Diameter_OC^3/6.0*d_oc)*10^6;
            %====================================================================================================
            result_core_oc_eabs=Mie(m_OC,Dx_wave_core);
            result_shell_OC=Mie(m_OC,Dx_wave_shell);
            Cabs_shell_OC(i)=result_shell_OC(3)*pi*power(Diameter_shell,2)/4/(pi*Diameter_shell^3/6.0*d_oc);
            Cabs_core_OC(i)= result_core_oc_eabs(3)*pi*power(Diameter_core(d),2)/4/(pi*Diameter_core(d)^3/6.0*d_oc);
            Eabsnonabs(i)=(1.8*BC_conc_factor*Cabs_core_BC_shell_OC(i)-1.2*shell_ratio^3*Cabs_shell_OC(i)+1.2*Cabs_core_OC(i))/(1.8*Cabs_core_BC(i)); 
            %============================================================================================================
            Cabs_shellratio_oc(i)=result_core_OC(3)*pi*power(Diameter_OC,2)/4/(pi*Diameter_OC^3/6.0*d_oc);
            Eabs_ocabs_ext(i)=(1.8*Cabs_core_BC(i)+1.2*(shell_ratio^3-1)*Cabs_shellratio_oc(i))/(1.8*Cabs_core_BC(i)); 
            
          end   
% ===========================================================纯雪  clear sky
data_out_shell=[];
for hh=1:numbers_hour(2);
    data_spec_shell_output=snicar8d_pw_stl(200,snow_radius(r),coszen(hh),0,1,data_bc_shell_spec_in,0,data_bc_OC_spec_in); 
    data_out_shell(:,1)=data_spec_shell_output(:,1);
    data_out_shell(:,hh+1)=data_spec_shell_output(:,2);
end
      %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@外混
%                 OC_conc=OC_conc_factor*BC_conc(n1);
for hh=1:numbers_hour(2);
                data_spec_core_output1=snicar8d_pw_stl(200,snow_radius(r),coszen(hh),bc_select(s),1,data_bc_core_spec_in,oc_select(s),data_bc_OC_spec_in);
                data_out_shell(:,numbers_hour(2)+1+hh)=data_spec_core_output1(:,2);
end
%===========================================================内混
for hh=1:numbers_hour(2);
                data_spec_shell_output1=snicar8d_pw_stl(200,snow_radius(r),coszen(hh),BC_conc_factor*bc_select(s),1,data_bc_shell_spec_in,0,data_bc_OC_spec_in);                 
                data_out_shell(:,2*numbers_hour(2)+1+hh)=data_spec_shell_output1(:,2);
%               data_spec_shell_averaged(:,n)=data_spec_shell_output(1:5,3);
end
data_out_shell_1=[];
data_out_shell_1=[data_out_shell_1 25 S_select S_select S_select];
data_out_shell_2=[];
data_out_shell_2=[data_out_shell_2 data_out_shell_1;data_out_shell];
            name=['_site_' num2str(site_select(s)) '_time_' num2str(R(ii)) '_snow_radius_' num2str(snow_radius(r))] ; % 
            eval(['save E:\2018ACP_pw\ACP\measure_data\clear_sky\BC_100\305_5000_spec_all_' name '.txt data_out_shell_2 -ascii;'] );%@@@@@@@@@@@@@@@@@
            %====================================
            % ===========================================================纯雪  cloudy sky
    data_out_shell=[];        
for hh=1:numbers_hour(2);
    data_spec_shell_output=snicar8d_pw_stl(200,snow_radius(r),coszen(hh),0,0,data_bc_shell_spec_in,0,data_bc_OC_spec_in); 
    data_out_shell(:,1)=data_spec_shell_output(:,1);
    data_out_shell(:,hh+1)=data_spec_shell_output(:,2);
end
      %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@外混
%                 OC_conc=OC_conc_factor*BC_conc(n1);
for hh=1:numbers_hour(2);
                data_spec_core_output1=snicar8d_pw_stl(200,snow_radius(r),coszen(hh),bc_select(s),0,data_bc_core_spec_in,oc_select(s),data_bc_OC_spec_in);
                data_out_shell(:,numbers_hour(2)+1+hh)=data_spec_core_output1(:,2);
end
%===========================================================内混
for hh=1:numbers_hour(2);
                data_spec_shell_output1=snicar8d_pw_stl(200,snow_radius(r),coszen(hh),BC_conc_factor*bc_select(s),0,data_bc_shell_spec_in,0,data_bc_OC_spec_in);                 
                data_out_shell(:,2*numbers_hour(2)+1+hh)=data_spec_shell_output1(:,2);
%               data_spec_shell_averaged(:,n)=data_spec_shell_output(1:5,3);
end
data_out_shell_1=[];
data_out_shell_1=[data_out_shell_1 25 S_select S_select S_select];
data_out_shell_2=[];
data_out_shell_2=[data_out_shell_2 data_out_shell_1;data_out_shell];
            name=['_site_' num2str(site_select(s)) '_time_' num2str(R(ii)) '_snow_radius_' num2str(snow_radius(r))] ; % 
            eval(['save E:\2018ACP_pw\ACP\measure_data\cloudy_sky\BC_100\305_5000_spec_all_' name '.txt data_out_shell_2 -ascii;'] );%@@@@@@@@@@@@@@@@@
            % =========================================================over
        end 
     end
   end
end
 disp('Hello World!');  
  %   
%             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\nonabs_shell\305_5000_spec_all_' name '.txt data_spec_shell_averaged -ascii;'] );

% =======                data_out_core1_averaged(:,n1)=data_spec_core_output1(1:5,3);
%     data_Eabs_otput(:,s+1)=Eabs_ocabs_ext;
%     data_Eabs_otput1(:,s+1)=Eabs;
%     data_Eabs_otput(:,1)=wave_length;
%     data_Eabs_otput1(:,1)=wave_length;
%             end;
%              name=['_Diameter_core_' num2str(Diameter_core(d)) '_snow_radius_' num2str(snow_radius(r))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
%             eval(['save E:\2018ACP_pw\ACP\measure_data\spectral_data\BC_OCshell\305_5000_spec_all_' name '.txt data_out_shell -ascii;'] )
%             %========================================================================================================
%             name1=['_Diameter_core_' num2str(Diameter_core(d)) '_snow_radius_' num2str(snow_radius(r))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
%             eval(['save E:\2018ACP_pw\ACP\measure_data\spectral_data\BC_OC\305_5000_spec_all_' name '.txt data_out_core1 -ascii;'] );
% %             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\BC_OCnonabs\305_5000_spec_all_' name '.txt data_out_core1_averaged -ascii;'] );
%         end
%              
%         %===================================================================================================================================== %   
%     %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++shell_ratio=0
%     %============================================================    
% 
%     name=['_Diameter_core_' num2str(Diameter_core(d))];  % '_Diameter_shell_' num2str(Diameter_shell)];
%     eval(['save E:\2018ACP_pw\ACP\measure_data\Eabs\BC_OC\305_5000_Eabs_all_' name '.txt data_Eabs_otput -ascii;']);
%     eval(['save E:\2018ACP_pw\ACP\measure_data\Eabs\BC_OCshell\305_5000_Eabs_all_' name '.txt data_Eabs_otput1 -ascii;']);
%     
% %     stop;
%    
% end;
disp('Hello World!');  
            
            