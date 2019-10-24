clc;
%eval('cd E:\2018ACP_pi\ACP\');
coszen=0.5;
bc_oc=csvread('E:\2018ACP_pw\ACP\OCDATA.csv',1,1);
oc_bc_ratio=bc_oc(:,1)./bc_oc(:,2);
index=find(oc_bc_ratio > 1 & oc_bc_ratio < 30);
number_sample=size(index,1);
site_select=bc_oc(index,3);
oc_bc_select=oc_bc_ratio(index);
bc_select=bc_oc(index,2);
oc_select=bc_oc(index,1);
bc_select=[0 ;bc_select];
oc_select=[0; oc_select];
oc_bc_select=[1; oc_bc_select];%252
%======================================
m_BC=complex(1.85,0.71);
% m_OC=complex(1.55,0.001);
d_bc=1.8;%  g cm-3
d_oc=1.2;
wave_length=305:10:4995;
k_OC_all=load('E:\2018ACP_pw\ACP\shitlwork\Mie_model\K_OC\300_750_K_OCk_brown_0.06_AAE_brown_-6_Diameter_cores_200nm.txt');
k_OC_all=[k_OC_all(6:10:451);zeros(1,425)'];
Diameter_core=[100 200 300 400 500]; %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@xiugai
% shell_ratio=[1.5 2 3 4];
size_wave=size(wave_length);
snow_radius=[100 200 500 1000];

    for d=1:5;
        for r=1:4;
            for s=1:number_sample+1;
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
%===========================================================ÄÚ»ì
                data_spec_shell_output=snicar8d_pw_stl(200,snow_radius(r),coszen,BC_conc_factor*bc_select(s),1,data_bc_shell_spec_in,0,data_bc_OC_spec_in);                 data_out_shell(:,1)=data_spec_shell_output(:,1);
                data_out_shell(:,1)=data_spec_shell_output(:,1);
                data_out_shell(:,s+1)=data_spec_shell_output(:,2);
%                 data_spec_shell_averaged(:,n)=data_spec_shell_output(1:5,3);
            
            name=['_Diameter_core_' num2str(Diameter_core(d)) '_snow_radius_' num2str(snow_radius(r))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
            eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result\nonabs_shell\305_5000_spec_all_' name '.txt data_out_shell -ascii;'] );%@@@@@@@@@@@@@@@@@
%             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\nonabs_shell\305_5000_spec_all_' name '.txt data_spec_shell_averaged -ascii;'] );%@@@@@@@@@@@
      %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Íâ»ì
%                 OC_conc=OC_conc_factor*BC_conc(n1);
                data_spec_core_output1=snicar8d_pw_stl(200,snow_radius(r),coszen,bc_select(s),1,data_bc_core_spec_in,oc_select(s),data_bc_OC_spec_in);
                data_out_core1(:,1)=data_spec_core_output1(:,1);
                data_out_core1(:,s+1)=data_spec_core_output1(:,2);
%                 data_out_core1_averaged(:,n1)=data_spec_core_output1(1:5,3);
    data_Eabs_otput(:,s+1)=Eabs_ocabs_ext;
    data_Eabs_otput1(:,s+1)=Eabs;
    data_Eabs_otput(:,1)=wave_length;
    data_Eabs_otput1(:,1)=wave_length;
            end;
             name=['_Diameter_core_' num2str(Diameter_core(d)) '_snow_radius_' num2str(snow_radius(r))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
            eval(['save E:\2018ACP_pw\ACP\measure_data\spectral_data\BC_OCshell\305_5000_spec_all_' name '.txt data_out_shell -ascii;'] )
            %========================================================================================================
            name1=['_Diameter_core_' num2str(Diameter_core(d)) '_snow_radius_' num2str(snow_radius(r))] ; % '_Diameter_shell_' num2str(Diameter_shell)];
            eval(['save E:\2018ACP_pw\ACP\measure_data\spectral_data\BC_OC\305_5000_spec_all_' name '.txt data_out_core1 -ascii;'] );
%             eval(['save E:\2018ACP_pw\ACP\SNICAR\snicar_result_spectral_averaged_0.3_1.0\BC_OCnonabs\305_5000_spec_all_' name '.txt data_out_core1_averaged -ascii;'] );
        end
             
        %===================================================================================================================================== %   
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++shell_ratio=0
    %============================================================    

    name=['_Diameter_core_' num2str(Diameter_core(d))];  % '_Diameter_shell_' num2str(Diameter_shell)];
    eval(['save E:\2018ACP_pw\ACP\measure_data\Eabs\BC_OC\305_5000_Eabs_all_' name '.txt data_Eabs_otput -ascii;']);
    eval(['save E:\2018ACP_pw\ACP\measure_data\Eabs\BC_OCshell\305_5000_Eabs_all_' name '.txt data_Eabs_otput1 -ascii;']);
    
%     stop;
   
end;
disp('Hello World!');  
            
            