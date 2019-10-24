clc;
   m_BC=complex(1.95,0.79);
   Diameter_core=10:10:1000;
   for i=1:100;
       Dx_wave_core=Diameter_core(i)*pi/550; %[qext qsca qabs qb asy qratio];
       result_core_BC=Mie(m_BC,Dx_wave_core);
       Cmee_core_BC(i)=result_core_BC(3)*pi*power(Diameter_core(i),2)/4/(pi*Diameter_core(i)^3/6.0*1.8)*10^3;%  6
   end
    p=semilogx(Diameter_core, Cmee_core_BC,'b');