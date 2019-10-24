function result=SNICAR_Eabs_stl_v(Diameter_core_in,Diameter_OC_in,shell_ratio_in,...
                                   snow_radius_in,BC_conc_in,abs_in)
% clc;
% clear;
coszen=0.5;
%data_spec_albedo=snicar8d_pw_stl(200,100,coszen,200,1);%密度 粒径 天顶角余弦值 BC浓度
%%直射1&漫射0   add(黑碳单次散射反照率 黑碳不对称因子 黑碳消光截面 ）利用米散射来求
m_BC=complex(1.95,0.79);
% m_OC=complex(1.55,10^-6);
d_bc=1.8;%  g cm-3
d_oc=1.2;
%data_bc_optic_in=zeros(470,3);
wave_length=305:10:4995;
%a=zeros(1,425);
abs=abs_in;
Diameter_core=Diameter_core_in;%100;%[100 200 300 400 500]; %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@xiugai
Diameter_OC=Diameter_OC_in;%200;
shell_ratio=shell_ratio_in;%[1.5 2 3];
k_nu=size(shell_ratio_in);
% shell_ratio=1.1:0.1:3;
size_wave=size(wave_length);
snow_radius=snow_radius_in;%[100 200 500 1000];
m_nu=size(snow_radius_in);
% snow_radius=50:50:1000;
BC_conc=BC_conc_in;%0:100:5000;
n_nu=size(BC_conc_in);
if Diameter_OC==100
    k_OC_all=load(['E:\2018ACP_pw\BioSNICAR_master\k_brown_0.3_in_550_AAE_6_Diameter_cores_100nm.txt']);
else
    k_OC_all=load(['E:\2018ACP_pw\BioSNICAR_master\k_brown_0.3_in_550_AAE_6_Diameter_cores_200nm.txt']);
end
k_OC_all=k_OC_all*0.5;
k_OC_all(k_OC_all<=10^-5)=10^-6; 
for j=1:1;
    for k=1:k_nu(2);
        Diameter_shell=shell_ratio(k)*Diameter_core(j);
        d_shell=(Diameter_core(j)/Diameter_shell)^3*d_bc+(1-(Diameter_core(j)/Diameter_shell)^3)*d_oc;
        BC_conc_factor=d_shell/d_bc*(Diameter_shell/Diameter_core(j))^3;
        OC_conc_factor=d_oc/d_bc*((Diameter_shell/Diameter_core(j))^3-1);
%         Diameter_OC=(Diameter_shell^3-Diameter_core(j)^3)^(1/3);
        for i=1:size_wave(2);
            Dx_wave_core=Diameter_core(j)*pi/wave_length(i);
            Dx_wave_shell=Diameter_shell*pi/wave_length(i);
            Dx_wave_OC=Diameter_OC*pi/wave_length(i);%@@@@@@@@@@@@@@@@@@@@@@@@
            if abs==1
                m_OC=complex(1.55,k_OC_all(i));
            else
                m_OC=complex(1.55,10^-6);
            end     
                      %@@@@@@@@@@@@@
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
        for m=1:m_nu(2);
            for n=1:n_nu(2);
                data_spec_shell_output=snicar8d_pw_stl_v(200,snow_radius(m),coszen,BC_conc_factor*BC_conc(n),1,data_bc_shell_spec_in,0,data_bc_OC_spec_in);                 data_out_shell(:,1)=data_spec_shell_output(:,1);
                data_out_shell(:,1)=data_spec_shell_output(:,1);
                data_out_shell(:,n+1)=data_spec_shell_output(:,2);
%                 data_spec_shell_averaged(:,n)=data_spec_shell_output(1:5,3);
            end;
            name=['_Diameter_core_' num2str(Diameter_core(j)) '_shell_ratio_' num2str(shell_ratio(k)) '_snow_radius_' num2str(snow_radius(m))] ; % '_Diameter_shell_' num2str(Diameter_shell)]
            if abs==1
                eval(['save E:\2018ACP_pw\ACP\SNICAR\data_omgea\OC_shell\305_5000_spec_all_' name '.txt data_out_shell -ascii;'] );%@@@@@@@@@@@@@@@@@
            else
                eval(['save E:\2018ACP_pw\ACP\SNICAR\data_omgea\nonabs_shell\305_5000_spec_all_' name '.txt data_out_shell -ascii;'] );
            end
                %             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\nonabs_shell\305_5000_spec_all_' name '.txt data_spec_shell_averaged -ascii;'] );%@@@@@@@@@@@
        end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@外混
    for m1=1:m_nu(2);
            for n1=1:n_nu(2);
                OC_conc=OC_conc_factor*BC_conc(n1);
                data_spec_core_output1=snicar8d_pw_stl_v(200,snow_radius(m1),coszen,BC_conc(n1),1,data_bc_core_spec_in,OC_conc,data_bc_OC_spec_in);
                data_out_core1(:,1)=data_spec_core_output1(:,1);
                data_out_core1(:,n1+1)=data_spec_core_output1(:,2);
%                 data_out_core1_averaged(:,n1)=data_spec_core_output1(1:5,3);
            end;
            name=['_Diameter_core_' num2str(Diameter_core(j)) '_shell_ratio_' num2str(shell_ratio(k)) '_snow_radius_' num2str(snow_radius(m1))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
            if abs==1
                eval(['save E:\2018ACP_pw\ACP\SNICAR\data_omgea\BC_OC\305_5000_spec_all_' name '.txt data_out_core1 -ascii;'] );
            else
                eval(['save E:\2018ACP_pw\ACP\SNICAR\data_omgea\BC_OCnonabs\305_5000_spec_all_' name '.txt data_out_core1 -ascii;'] );
            end
            omgeaaaa=m1
%             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\BC_OCnonabs\305_5000_spec_all_' name '.txt data_out_core1_averaged -ascii;'] );
    end
         
%         data_Eabs_otput(:,k+1)=Eabs_ocabs_ext;
        %===================================================================================================================================== %   
 
%======================================================================================================================================
    end;    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++shell_ratio=0
       for m2=1:m_nu(2);
            for n2=1:n_nu(2);
                data_spec_core_output=snicar8d_pw_stl_v(200,snow_radius(m2),coszen,BC_conc(n2),1,data_bc_core_spec_in,0,data_bc_OC_spec_in);
                data_out_core(:,1)=data_spec_core_output(:,1);
                data_out_core(:,n2+1)=data_spec_core_output(:,2);
%                 data_out_core_averaged(:,n2)=data_spec_core_output(1:5,3);
            end;
            name=['_Diameter_core_' num2str(Diameter_core(j)) '_shell_ratio_' num2str(0) '_snow_radius_' num2str(snow_radius(m2))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
            eval(['save E:\2018ACP_pw\ACP\SNICAR\data_omgea\BC_OC\305_5000_spec_all_' name '.txt data_out_core -ascii;'] ); %@@@@@@@@@@@@@@@@@@
			eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar1.5\BC_OMEGA\305_5000_spec_all_' name '.txt data_out_core -ascii;'] );
%             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\BC_OCnonabs\305_5000_spec_all_' name '.txt data_out_core_averaged -ascii;'] );
       end;
	   
	   %=======================================================
	   for m2=1:m_nu(2);
            for n2=1:n_nu(2);
                data_spec_core_output=snicar8d_pw_stl_v(200,snow_radius(m2),coszen,1.5*BC_conc(n2),1,data_bc_core_spec_in,0,data_bc_OC_spec_in);
                data_out_core(:,1)=data_spec_core_output(:,1);
                data_out_core(:,n2+1)=data_spec_core_output(:,2);
%                 data_out_core_averaged(:,n2)=data_spec_core_output(1:5,3);
            end;
            name=['_Diameter_core_' num2str(Diameter_core(j)) '_shell_ratio_' num2str(0) '_snow_radius_' num2str(snow_radius(m2))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
%           eval(['save E:\2018ACP_pw\ACP\SNICAR\data_omgea\BC_OCnonabs\305_5000_spec_all_' name '.txt data_out_core -ascii;'] ); %@@@@@@@@@@@@@@@@@@
			eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar1.5\BC_OMEGA1.5\305_5000_spec_all_' name '.txt data_out_core -ascii;'] );
%             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\BC_OCnonabs\305_5000_spec_all_' name '.txt data_out_core_averaged -ascii;'] );
       end;

    %============================================================    
%     data_Eabs_otput(:,1)=wave_length;
%     name=['_Diameter_core_' num2str(Diameter_core(j))];  % '_Diameter_shell_' num2str(Diameter_shell)];
%     eval(['save E:\2018ACP_pw\ACP\SNICAR\mie_result\BC_OC\305_5000_Eabs_all_' name '.txt data_Eabs_otput -ascii;']);
%     stop;
   
end;
result='Hello World!';
%disp('Hello World!');  


% 
%    m_BC=complex(1.85,0.71);
%    Diameter_core=10:10:1000;
%    for i=1:100;
%        Dx_wave_core=Diameter_core(i)*pi/550; %[qext qsca qabs qb asy qratio];
%        result_core_BC=Mie(m_BC,Dx_wave_core);
%        Cmee_core_BC(i)=result_core_BC(3)/result_core_BC(1)*pi*power(Diameter_core(i),2)/4/(pi*Diameter_core(i)^3/6.0*1.8)*10^3;%  6
%    end
%     p=semilogx(Diameter_core, Cmee_core_BC,'b');
%    result_core_BC;
%    result_core_BC_shell_OC;
% %data_bc_optic_in(:,1)=wave_legth;
% %data_spec_albedo=snicar8d_pw_stl(200,100,coszen,200,1,data_bc_optic_in);
% 
% clear;   
% % result=SizeDist_Optics(1.55+0.71i, 200, 0.1, 550);
% % result;
% m_BC=complex(1.85,0.71); %550nm 380nm
% % m_BC=complex(2.0,0.71); %550nm 380nm
% %m_OC=complex(1.55,0.06); %380nm
% % wavelength=300:10:750;
% % wavelength=[365 532];
% % wavelength_all=365:750;
% wavelength=350:750;
% % index_532=find(wavelength_all==532);
% % a=wavelength_all(index_532(1))
% % stop
% size_wave=size(wavelength);
% m_nonabs=complex(1.55,0);
% k_brown=[0.03 0.06 0.09];
% AAE_brown=[-4,-6];
% Diameter_core=100;
% Diameter_shell=200;
% k_OC_all=load(['E:\2018ACP_pw\ACP\shitlwork\Mie_model\K_OC\350_750_K_OCk_brown_0.06_AAE_brown_-6_Diameter_cores_200nm.txt']);
% 
% % stop
% for i=1:size_wave(2);
%     Dx_wave_core=Diameter_core*pi/wavelength(i);
%     Dx_wave_shell=Diameter_shell*pi/wavelength(i);
%     m_OC=complex(1.55,k_OC_all(i));
%     result_core = Mie(m_OC,Dx_wave_core);
%     result_shell = Mie(m_OC,Dx_wave_shell);
%     Cabs_OC_difference(i)=result_shell(3)*pi*power(Diameter_shell,2)/4-result_core(3)*pi*power(Diameter_core,2)/4;
%     result_core_BC= Mie(m_BC,Dx_wave_core);
%     Cabs_core_BC(i)=result_core_BC(3)*pi*power(Diameter_core,2)/4;
%     result_core_BC_shell_OC=Miecoated(m_BC,m_OC,Dx_wave_core,Dx_wave_shell,1);
%     Cabs_core_BC_shell_OC(i)=result_core_BC_shell_OC(3)*pi*power(Diameter_shell,2)/4;
%     result_core_BC_shell_nonabs=Miecoated(m_BC,m_nonabs,Dx_wave_core,Dx_wave_shell,1);
%     Cabs_core_BC_shell_nonabs(i)=result_core_BC_shell_nonabs(3)*pi*power(Diameter_shell,2)/4;
% %     stop;
%    
% end
% 
% Cabs_core_BC_lensing=Cabs_core_BC_shell_OC-Cabs_OC_difference;
% Eabs=Cabs_core_BC_lensing./Cabs_core_BC;  
% Eabs_nonabs=Cabs_core_BC_shell_nonabs./Cabs_core_BC;
% Eabs_lense_OC=Cabs_core_BC_shell_OC./Cabs_core_BC;
% Eabs_nolense=(Cabs_OC_difference+Cabs_core_BC)./Cabs_core_BC;
% p=plot(wavelength,Eabs_nonabs,'k');
% hold on;
% p=plot(wavelength,Eabs,'r');
% hold on;
% % for i=1:11
% p=plot(wavelength,Eabs_lense_OC,'b');
% hold on;
% % end
% set(gca,'xlim',[350 600],'ylim',[0 3]);
% Eabs_all(:,1)=wavelength;
% Eabs_all(:,2)=Eabs;
% Eabs_all(:,3)=Eabs_lense_OC;
% Eabs_all(:,4)=Eabs_nonabs;
% Eabs_all(:,5)=Eabs_nolense;
% name=['_Diameter_core_' num2str(Diameter_core) '_Diameter_shell_' num2str(Diameter_shell)];
% eval(['save  E:\2018ACP_pw\ACP\shitlwork\Mie_model\K_OC\Eabs\350_750_Eabs_all_' name '.txt Eabs_all -ascii;'] );
% %save F:\work\2017\ACP\Mie_model\K_OC\Eabs\Eabs_all_K_OCk_brown_0.06_AAE_brown_-6_Diameter_cores_200nm.txt Eabs_all -ascii;
% 
% Diameter=200; %nm
% wavelength=380; %nm
% Density_BC=1.8; %g/cm3
% Dx_wave=Diameter*pi/wavelength;   %  ;Particle size parameter(s)
% for i=1:1500;
%     m_OC=complex(1.55,i*0.0001);
%     Dx=Dx_wave;
%     result = Mie(m_OC,Dx);
%     ssa(i)=result(2)/result(1);
%     k_OC(i)=i*0.0001;
% %     MAC(i)=result(3)*3.0/4.0/((Diameter(i)/2.0)*1.8)*1000;
% end
% ssa=abs(ssa-0.75);
% index=find(ssa==min(ssa));
% k_OC(index)
% p=plot(Diameter,MAC,'xlim',[10 1000],'ylim',[0 7],'xscale','log');
% plot(Diameter,MAC);
% set(gca,'xlim',[10 1000],'ylim',[0 7],'xscale','log');
