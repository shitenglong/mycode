clear;
BC_con=[0:10:100 120:20:500 550:50:1000 1100:100:2000 2200:200:5000 5500:500:10000];% 10000:1000:50000];
size_BC=size(BC_con);
Solar_zenith=20:0.5:75;
size_SZ=size(Solar_zenith);
Snow_size=[10:5:200 210:10:500 520:20:1000];
size_snowsize=size(Snow_size);
Cos_SZ=cos(Solar_zenith/180*pi);
% stop
load 'F:\gr010318_000.txt';
load 'F:\gr010318_004.txt';
data_out = snicar8d_pw(0.005,300,170,cos(60/180*pi),100,0);
hold on;
p=plot(gr010318_004(:,1)/1000,gr010318_004(:,4)/100,'r');
hold on;
%  data_out = snicar8d_pw(0.5,300,Snow_size(k),Cos_SZ(j),BC_con(i))
% stop
% for  i=1: size_BC(2)
%     ss=i
%     number=1;
%     for j=1:size_SZ(2);
%         j=11;
%         xx=j
%         for k=1:size_snowsize(2);
%             k=14;
%             data_out = snicar8d_pw(0.5,300,Snow_size(k),Cos_SZ(j),BC_con(i)); %snow_thickness (m), snow_density (kg/m3), snow_size_radius (um),cos(Solar_zenith),BC_conc (ng g-1)
%             hold on;
%             data_out_albedo(1:220,number)=data_out(1:220,2);
%             data_wave=data_out(1:220,1);
%             number=number+1;
%               name=['Snicar_wave'];
%                 eval(['save F:\work\model\snow_albedo_model\SNICAR\BioSNICAR_master\BioSNICAR_master\SNICAR_out\Albedo_' name '.txt data_wave -ascii;'] );   
%             stop
%         end;
%         stop
%     end
%     stop
%     name=['BC_con_' num2str(BC_con(i))];
%     eval(['save F:\work\model\snow_albedo_model\SNICAR\BioSNICAR_master\BioSNICAR_master\SNICAR_out\Albedo_' name '.txt data_out_albedo -ascii;'] );   

% end

% index=find(data_out(:,1)==1.245);
% BC_i=BC_con(i)
% Albedo_1245=data_out(index,2)
% hold on;
% stop
% end;
%   p=plot(Snicar_online_7(:,1),Snicar_online_7(:,2),'r');
% hold on;