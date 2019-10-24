clc;
clear;
coszen=0.5;
%data_spec_albedo=snicar8d_pw_stl(200,100,coszen,200,1);%密度 粒径 天顶角余弦�? BC浓度
%%直射1&漫射0   add(黑碳单次散射反照�?黑碳不对称因�?黑碳消光截面 ）利用米散射来求
m_BC=complex(1.95,0.79);
m_OC=complex(1.55,10^-6);
d_bc=1.8;%  g cm-3
d_oc=1.2;
%data_bc_optic_in=zeros(470,3);
wave_length=305:10:4995;
%a=zeros(1,425);
k_OC_all=load(['E:\2018ACP_pw\BioSNICAR_master\k_brown_0.3_in_550_AAE_6_Diameter_cores_200nm.txt']);
k_OC_all(k_OC_all<=10^-5)=10^-6; 
% k_OC_all=[k_OC_all(6:10:451);zeros(1,425)'];
Diameter_core=200; %[100 200 300 400 500]; %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@xiugai
Diameter_OC=200;                           %@@@@@@@@@@@@@@@@@@@@@@@@@XIUGAI
 shell_ratio=[1.5 2 2.5];
%shell_ratio=1.1:0.1:3;
size_wave=size(wave_length);
snow_radius=[100 200 300 500];        %%@@@@@@@@@@@@@@@@@@@@
%snow_radius=50:50:1000;
BC_conc=0:100:5000;
for j=1:1;
    for k=1:3;
        Diameter_shell=shell_ratio(k)*Diameter_core(j);
        d_shell=(Diameter_core(j)/Diameter_shell)^3*d_bc+(1-(Diameter_core(j)/Diameter_shell)^3)*d_oc;
        BC_conc_factor=d_shell/d_bc*(Diameter_shell/Diameter_core(j))^3;
        OC_conc_factor=d_oc/d_bc*((Diameter_shell/Diameter_core(j))^3-1); %@@@@@@@@@@@@@@@@@error
%         Diameter_OC=(Diameter_shell^3-Diameter_core(j)^3)^(1/3);
        for i=1:size_wave(2);
            Dx_wave_core=Diameter_core(j)*pi/wave_length(i);
            Dx_wave_shell=Diameter_shell*pi/wave_length(i);
            Dx_wave_OC=Diameter_OC*pi/wave_length(i);%@@@@@@@@@@@@@@@@@@@@@@@@
%             m_OC=complex(1.55,k_OC_all(i));          %@@@@@@@@@@@@@
           % result_core = Mie(m_OC,Dx_wave_core);
           % result_shell = Mie(m_OC,Dx_wave_shell);
           % Cabs_OC_difference(i)=result_shell(3)*pi*power(Diameter_shell,2)/4-result_core(3)*pi*power(Diameter_core,2)/4;
            result_core_BC=Mie(m_BC,Dx_wave_core);
            Cabs_core_BC(i)=result_core_BC(3)*pi*power(Diameter_core(j),2)/4/(pi*Diameter_core(j)^3/6.0*d_bc);
            result_core_BC_shell_OC=Miecoated(m_BC,m_OC,Dx_wave_core,Dx_wave_shell,1);
            Cabs_core_BC_shell_OC(i)=result_core_BC_shell_OC(3)*pi*power(Diameter_shell,2)/4/(pi*Diameter_shell^3/6.0*d_shell);
            Eabs(i)=BC_conc_factor*Cabs_core_BC_shell_OC(i)/Cabs_core_BC(i); 
            data_bc_core_spec_in(i,1)=result_core_BC(2)/result_core_BC(1);
            data_bc_core_spec_in(i,2)=result_core_BC(5);
            data_bc_core_spec_in(i,3)=result_core_BC(1)*pi*power(Diameter_core(j),2)/4/(pi*Diameter_core(j)^3/6.0*d_bc)*10^6;
            %result_core_BC_shell_nonabs=Miecoated(m_BC,m_nonabs,Dx_wave_core,Dx_wave_shell,1);
            %Cabs_core_BC_shell_nonabs(i)=result_core_BC_shell_nonabs(3)*pi*power(Diameter_shell,2)/4;
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
            Cabs_core_OC(i)= result_core_oc_eabs(3)*pi*power(Diameter_core(j),2)/4/(pi*Diameter_core(j)^3/6.0*d_oc);
            Eabsnonabs(i)=(1.8*BC_conc_factor*Cabs_core_BC_shell_OC(i)-1.2*shell_ratio(k)^3*Cabs_shell_OC(i)+1.2*Cabs_core_OC(i))/(1.8*Cabs_core_BC(i)); 
            %============================================================================================================
            Cabs_shellratio_oc(i)=result_core_OC(3)*pi*power(Diameter_OC,2)/4/(pi*Diameter_OC^3/6.0*d_oc);
            Eabs_ocabs_ext(i)=(1.8*Cabs_core_BC(i)+1.2*(shell_ratio(k)^3-1)*Cabs_shellratio_oc(i))/(1.8*Cabs_core_BC(i)); 
            
        end   
%===========================================================内混
%         for m=1:4;
%             for n=1:51;
%                 data_spec_shell_output=snicar8d_pw_stl(200,snow_radius(m),coszen,BC_conc_factor*BC_conc(n),1,data_bc_shell_spec_in,0,data_bc_OC_spec_in);                 data_out_shell(:,1)=data_spec_shell_output(:,1);
%                 data_out_shell(:,1)=data_spec_shell_output(:,1);
%                 data_out_shell(:,n+1)=data_spec_shell_output(:,2);
% %                 data_spec_shell_averaged(:,n)=data_spec_shell_output(1:5,3);
%             end;
%             name=['_Diameter_core_' num2str(Diameter_core(j)) '_shell_ratio_' num2str(shell_ratio(k)) '_snow_radius_' num2str(snow_radius(m))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
%             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_v1\OC_shell\305_5000_spec_all_' name '.txt data_out_shell -ascii;'] );%@@@@@@@@@@@@@@@@@
% %             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\nonabs_shell\305_5000_spec_all_' name '.txt data_spec_shell_averaged -ascii;'] );%@@@@@@@@@@@
%         end
%    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@外混
%     for m1=1:20;
%             for n1=1:1001;
%                 OC_conc=OC_conc_factor*BC_conc(n1);
%                 data_spec_core_output1=snicar8d_pw_stl(200,snow_radius(m1),coszen,BC_conc(n1),1,data_bc_core_spec_in,OC_conc,data_bc_OC_spec_in);
%                 data_out_core1(:,1)=data_spec_core_output1(:,1);
%                 data_out_core1(:,n1+1)=data_spec_core_output1(:,2);
% %                 data_out_core1_averaged(:,n1)=data_spec_core_output1(1:5,3);
%             end;
%             name=['_Diameter_core_' num2str(Diameter_core(j)) '_shell_ratio_' num2str(shell_ratio(k)) '_snow_radius_' num2str(snow_radius(m1))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
%             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_v1\BC_OC\305_5000_spec_all_' name '.txt data_out_core1 -ascii;'] );
% %             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\BC_OCnonabs\305_5000_spec_all_' name '.txt data_out_core1_averaged -ascii;'] );
%     aaa_abs=m1
%     end
         
        data_Eabs_otput(:,k+1)=Eabs./Eabs_ocabs_ext;
        %===================================================================================================================================== %   
    for m2=1:4;
            for n2=1:51;
                data_spec_core_output=snicar8d_pw_stl(200,snow_radius(m2),coszen,1.5*BC_conc(n2),1,data_bc_core_spec_in,0,data_bc_OC_spec_in);
                data_out_core(:,1)=data_spec_core_output(:,1);
                data_out_core(:,n2+1)=data_spec_core_output(:,2);
%                 data_out_core_averaged(:,n2)=data_spec_core_output(1:5,3);
            end;
            name=['_Diameter_core_' num2str(Diameter_core(j)) '_shell_ratio_' num2str(0) '_snow_radius_' num2str(snow_radius(m2))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
            eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar1.5\BC1.5\305_5000_spec_all_' name '.txt data_out_core -ascii;'] ); %@@@@@@@@@@@@@@@@@@
%             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\BC_OCnonabs\305_5000_spec_all_' name '.txt data_out_core_averaged -ascii;'] );
    end;
%======================================================================================================================================
figure;
plot(wave_length,Cabs_core_BC_shell_OC*10^3)
xlim([300 1400])
title(['BC\_100nm\_nonabsshell\_ratio\_' num2str(shell_ratio(k))])
xlabel('wavelength')
ylabel('MAC (g m^{-2})')
saveas(gcf,['E:\2018ACP_pw\MAC_FIGURE\Diameter_BC_' num2str(Diameter_core) '_absshell_ratio_' num2str(shell_ratio(k)) '.png']);
% saveas(gcf,['E:\2018ACP_pw\MAC_FIGURE\Diameter_OC_' num2str(Diameter_OC) '.png']);
close

    end;    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++shell_ratio=0

    %============================================================ 吸收增强   
    data_Eabs_otput(:,1)=wave_length;
    name=['_Diameter_core_' num2str(Diameter_core(j))];  % '_Diameter_shell_' num2str(Diameter_shell)];
    eval(['save E:\2018ACP_pw\ACP\SNICAR\mie_result\nonabs_shell\305_5000_Eabs_all_' name '.txt data_Eabs_otput -ascii;']);
%     stop;
   
end;
disp('Hello World!');




