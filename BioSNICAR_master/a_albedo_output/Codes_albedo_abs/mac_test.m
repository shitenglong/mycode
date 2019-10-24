clc;
clear;
wave_length=305:10:5000;%380:750;%
k_OC_all=load(['E:\2018ACP_pw\BioSNICAR_master\k_brown_0.3_in_550_AAE_6_Diameter_cores_200nm.txt']);
k_OC_all(k_OC_all<=10^-5)=10^-6; 
Diameter_OC=200;
d_oc=1.2;
m_BC=complex(1.95,0.79);
Diameter_BC=150;
d_bc=1.8;
size_wave=size(wave_length);

for i=1:size_wave(2);
            Dx_wave_OC=Diameter_OC*pi/wave_length(i);%@@@@@@@@@@@@@@@@@@@@@@@@
            m_OC=complex(1.55,k_OC_all(i)); 
            Dx_wave_bc=Diameter_BC*pi/wave_length(i);
            result_core_BC=Mie(m_BC,Dx_wave_bc);
            result_core_OC=Mie(m_OC,Dx_wave_OC);
            Cabs_core_BC(i)=result_core_BC(3)*pi*power(Diameter_BC,2)/4/(pi*Diameter_BC^3/6.0*d_bc)*10^3;
            Cabs_core_OC(i)=result_core_OC(3)*pi*power(Diameter_OC,2)/4/(pi*Diameter_OC^3/6.0*d_oc)*10^3;
end
figure;
plot(wave_length,Cabs_core_BC)
xlim([300 1400])
title('BC 150nm')
xlabel('wavelength')
ylabel('MAC (g m^{-2})')

% saveas(gcf,['E:\2018ACP_pw\MAC_FIGURE\Diameter_BC_' num2str(Diameter_BC) '_shell_ratio_' num2str(0) '.png']);
% saveas(gcf,['E:\2018ACP_pw\MAC_FIGURE\Diameter_OC_' num2str(Diameter_OC) '.png']);
% close


% wavelength=305:10:5000;%380:750;%
% size_wave=size(wavelength);
% abs_550_oc=0.3  ;  %0.002
% abs_oc=exp(-6.0*log(wavelength./550.0))*abs_550_oc;
%          Diameter_OC=200;
%          d_oc=1.2;
%          K_brown_range=0.00001:0.0001:0.5;
%          size_wave_k_range=size(K_brown_range);
%         for i=1:size_wave(2);
%             for j=1:size_wave_k_range(2);
%             Dx_wave_core=Diameter_OC*pi/wavelength(i);
%             m_oc_380=complex(1.55,K_brown_range(j));
%             result=Mie(m_oc_380, Dx_wave_core);
%             qabs_j(j)=result(3)*pi*power(Diameter_OC,2)/4/(pi*Diameter_OC^3/6.0*d_oc)*10^3;
%             end
%             qabs_j_bias=abs(qabs_j-abs_oc(i));
%             index=find(qabs_j_bias==min(qabs_j_bias));
%             K_oc(i)=K_brown_range(index(1));
% %             stop
%         end
%         p=plot(wavelength,K_oc,'b');
%         set(gca,'xlim',[300,750],'ylim',[0 0.4]);   
%         hold on;
%         p=plot(wavelength,K_oc,'r');
%         K_oc_1=K_oc'
%         name=['k_brown_0.3_in_550_AAE_6_Diameter_cores_200nm'];
%         eval(['save E:\2018ACP_pw\' name '.txt K_oc_1 -ascii;'] );



%=============================================================bc fit
clc;
clear;
wave_length=405:10:750;%380:750;%
m_BC=complex(1.95,0.79);
Diameter_BC=150;
d_bc=1.8;
size_wave=size(wave_length);
for i=1:size_wave(2);
            Dx_wave_bc=Diameter_BC*pi/wave_length(i);
            result_core_BC=Mie(m_BC,Dx_wave_bc);
            Cabs_core_BC(i)=result_core_BC(3)*pi*power(Diameter_BC,2)/4/(pi*Diameter_BC^3/6.0*d_bc)*10^3;
end
xx_fit=log(wave_length/455.0);
yy_fit=log(Cabs_core_BC/8.5320);
fun=@(a,x)a*x; 
a=lsqcurvefit(fun,0,xx_fit,yy_fit);%»Ø¹éÖ±Ïßx1=0:0.1:10;y1=a*x1;hold on; plot(x1,y1,'r-');

figure;
plot(wave_length,Cabs_core_BC)
xlim([300 1400])
title('BC 150nm')
xlabel('wavelength')
ylabel('MAC (g m^{-2})')
hold on;
xx=305:10:5000;
yy=(xx/455.0).^a*8.5320;
plot(xx,yy);
hold on
         r_bc_range=0.1:0.01:2;
         size_r_bc_range=size(r_bc_range);
         K_brown_range=0.1:0.01:2;
         size_wave_k_range=size(K_brown_range);
        for i=1:15;
            for k=1:size_r_bc_range(2);
            for j=1:size_wave_k_range(2);
            Dx_wave_core=Diameter_BC*pi/xx(i);
            m_oc_380=complex(r_bc_range(k),K_brown_range(j));
            result=Mie(m_oc_380, Dx_wave_core);
            qabs_j(k,j)=result(3)*pi*power(Diameter_BC,2)/4/(pi*Diameter_BC^3/6.0*d_bc)*10^3;
            end
            end
            qabs_j_bias=abs(qabs_j-yy(i));
            aaa=min(qabs_j_bias);
%             bbb=min(aaa);
           [row,column]=find(qabs_j_bias==min(aaa));
            r_oc(i)=r_bc_range(row);
            K_oc(i)=K_brown_range(column);
%             stop
        end
        for i=11:470
             K_oc(i)=0.71;
        end
        
        for i=1:470;
            Dx_wave_bc1=Diameter_BC*pi/xx(i);
            m_BC1=complex(1.95,K_oc(i));
            result_core_BC1=Mie(m_BC1,Dx_wave_bc1);
            Cabs_core_BC1(i)=result_core_BC1(3)*pi*power(Diameter_BC,2)/4/(pi*Diameter_BC^3/6.0*d_bc)*10^3;
        end
        plot(xx,Cabs_core_BC1)

%==============================================================


%    m_BC=complex(1.95,0.79);
%    Diameter_core=10:10:1000;
%    for i=1:100;
%        Dx_wave_core=Diameter_core(i)*pi/550; %[qext qsca qabs qb asy qratio];
%        result_core_BC=Mie(m_BC,Dx_wave_core);
%        Cmee_core_BC(i)=result_core_BC(3)*pi*power(Diameter_core(i),2)/4/(pi*Diameter_core(i)^3/6.0*1.8)*10^3;%  6
%    end
%     p=semilogx(Diameter_core, Cmee_core_BC,'b');