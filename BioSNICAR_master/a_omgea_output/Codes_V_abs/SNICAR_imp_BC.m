clear;
BC_con=0:5:2000;% 10000:1000:50000];
size_BC=size(BC_con);
Solar_zenith=[50 60 70];
size_SZ=size(Solar_zenith);
Snow_size=100:100:500;
size_snowsize=size(Snow_size);
Cos_SZ=cos(Solar_zenith/180*pi);
for i=1:size_snowsize(2);
    for j=1:size_SZ(2);
        for k=1:size_BC(2);
%             k=k
            data_out = snicar8d_pw(0.5,300,Snow_size(i),Cos_SZ(j),BC_con(k));%snow_thickness (m), snow_density (kg/m3), snow_size_radius (um),cos(Solar_zenith),BC_conc (ng g-1)
            
            imp(k)=log(data_out(25,2))/log(data_out(94,2));
            Albedo(:,k)=data_out(1:220,2);
            %             wave=data_out(:,1);
%             index_550=find(wave==0.5550);
%             index_1240=find(wave==1.235);
%             stop
        end
%         p=plot(BC_con,imp);
%         hold on;
%             Albedo1=Albedo';
            name=['Oradius_' num2str(Snow_size(i)) '_SZA_' num2str(Solar_zenith(j))];
            eval(['save F:\work\2017\RSE\Northeastern_China\MYD09GA\Figure1\SNICAR\snow_albedo_reduction_BC\snow_albedo_reduction_BC_' name '.txt Albedo -ascii;'] );
%             stop
    end
end