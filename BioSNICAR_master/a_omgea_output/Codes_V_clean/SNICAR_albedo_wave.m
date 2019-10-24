clear;
% BC_con=0:5:2000;% 10000:1000:50000];
BC_con=[0 100 1000];
size_BC=size(BC_con);
Solar_zenith=60;
size_SZ=size(Solar_zenith);
Snow_size=100:100:300;
size_snowsize=size(Snow_size);
Cos_SZ=cos(Solar_zenith/180*pi);
% stop
 for k=1:size_BC(2);
    for j=1:size_SZ(2);
       for i=1:size_snowsize(2);
%             k=k
            data_out = snicar8d_pw(0.5,300,Snow_size(i),Cos_SZ(j),BC_con(k));%snow_thickness (m), snow_density (kg/m3), snow_size_radius (um),cos(Solar_zenith),BC_conc (ng g-1)
            
            imp(k)=log(data_out(25,2))/log(data_out(94,2));
            Albedo(:,i)=data_out(1:220,2);
            %             wave=data_out(:,1);
%             index_550=find(wave==0.5550);
%             index_1240=find(wave==1.235);
%             stop
        end
%         p=plot(BC_con,imp);
%         hold on;
%             Albedo1=Albedo';
            name=['BC_con_' num2str(BC_con(k))];
            eval(['save F:\work\2017\RSE\Northeastern_China\MYD09GA\Figure1\SNICAR\snow_albedo_wavelength\snow_albedo_wavelengths_' name '.txt Albedo -ascii;'] );
%             stop
    end
end