% Snow_size=[10:5:200 210:10:500 520:20:1000];
Snow_size=[10:5:200 210:10:500];
size_snowsize=size(Snow_size);
load 'F:\work\2015ACP\reflect\BCbest_Fe2.8%_meansnowdepth\Zhangmodel_MAC_single_layer\measured_albedo\divided_by_1.19\Measured_albedos_1240nm.txt'

for i=1:6;
    BC_con=0;
    SZA=Measured_albedos_1240nm(i,3);
    Cos_SZ=cos(SZA/180*pi);
    for j=1:2;
        for k=1:size_snowsize(2);
%             k=14;
            data_out = snicar8d_pw(0.5,300,Snow_size(k),Cos_SZ,BC_con); %snow_thickness (m), snow_density (kg/m3), snow_size_radius (um),cos(Solar_zenith),BC_conc (ng g-1)
            Albedo_1240(k,1:2)=data_out(94,1:2);
            
            %             hold on;
%             stop
%             data_out_albedo(1:220,number)=data_out(1:220,2);
        end
      Albedo_1240_abs=abs(Albedo_1240(:,2)-Measured_albedos_1240nm(i,j));
      index=find(Albedo_1240_abs==min(Albedo_1240_abs));
      
      R_snicar(i,j)=Snow_size(index(1));
    end
end