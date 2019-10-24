function HeightAngle = SolarAngle(year_in,month_in,day_in,lon_in,lat_in)
% function [HeightAngle AzimuthAngle] = SolarAngle(DateTime,Loaction)
% 功能：计算太阳天顶距和太阳方位角 SolarAngle.cpp
% 作者：汪帮主 2009-6-6 wzj08@mails.jlu.edu.cn 	
% 修改：liuzhiping 2010-12-3 zhpliu@cumt.edu.cn
% input DateTime = [year month day hour min sec];
% input Loaction = [longitude latitude TimeZone];

year = year_in;
month = month_in;
day = day_in;
hour = 1:1:23;
min = 0;
sec = 0;
lon = lon_in;
lat = lat_in;
TimeZone = 0;

% //年积日的计算
JD1 = fix(365.25*(year-1))+fix(30.6001*(1+13))+1+hour/24+1720981.5;
if month<=2
     JD2 = fix(365.25*(year-1))+fix(30.6001*(month+13))+day+hour/24+1720981.5;
else
    JD2 = fix(365.25*year)+fix(30.6001*(month+1))+day+hour/24+1720981.5;
end
DOY = JD2-JD1+1;
N0 = 79.6764 + 0.2422*(year-1985) - floor((year-1985)/4.0);
sitar = 2*pi*(DOY-N0)/365.2422;
ED = 0.3723 + 23.2567*sin(sitar) + 0.1149*sin(2*sitar) - 0.1712*sin(3*sitar)- 0.758*cos(sitar) + 0.3656*cos(2*sitar) + 0.0201*cos(3*sitar);
ED = ED*pi/180;           % //ED本身有符号，无需判断正负。

if (lon >= 0)	
    if (TimeZone == -13)
        dLon = lon - (floor((lon*10-75)/150)+1)*15.0;
    else
        dLon = lon - TimeZone*15.0;   % //地球上某一点与其所在时区中心的经度差
    end
else 
    if (TimeZone == -13)
        dLon =  (floor((lon*10-75)/150)+1)*15.0- lon;
    else
        dLon =  TimeZone*15.0- lon;
    end
end
dLon=0
Et = 0.0028 - 1.9857*sin(sitar) + 9.9059*sin(2*sitar) - 7.0924*cos(sitar)- 0.6882*cos(2*sitar); % //视差

gtdt = hour + min/60.0 + sec/3600.0 + dLon/15; %//地方时
gtdt = gtdt + Et/60.0;
dTimeAngle = 15.0*(gtdt-12);
dTimeAngle = dTimeAngle*pi/180;	
latitudeArc = lat*pi/180;

HeightAngleArc = asin(sin(latitudeArc)*sin(ED)+cos(latitudeArc)*cos(ED)*cos(dTimeAngle)); % //高度角计算公式

CosAzimuthAngle = (sin(HeightAngleArc)*sin(latitudeArc)-sin(ED))/cos(HeightAngleArc)/cos(latitudeArc); % //方位角计算公式
AzimuthAngleArc = acos(CosAzimuthAngle);
HeightAngle = HeightAngleArc*180/pi;
AzimuthAngle = AzimuthAngleArc *180/pi;

if( dTimeAngle < 0 )
    AzimuthAngle = 180 - AzimuthAngle;
else
    AzimuthAngle = 180 + AzimuthAngle;
end


