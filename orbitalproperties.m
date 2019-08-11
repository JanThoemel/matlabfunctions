function [density,v,meanRadiusOfEarth,mu,MeanMotion]=orbitalproperties(altitude)

  meanRadiusOfEarth=6371000;      %% [m]
  equatorialRadiusOfEarth=6378.1; %% [m]
  polarRadiusOfEarth=6356.8;      %% [m] 
  mu=3.986004418E14;              %% [m3?s?2] gravitational constant

  densityarray=  [300.0 2.458E-14;
                  350.0 9.025E-15;
                  400.0 3.560E-15;
                  450.0 1.472E-15;
                  500.0 6.304E-16;
                  550.0 2.784E-16;
                  600.0 1.270E-16];
  densityarray(:,1)=densityarray(:,1)*1000;         %from km to m
  densityarray(:,2)=densityarray(:,2)/1000*100^3;   %from g/cm3 to kg/m3
  f=fit(densityarray(:,1),densityarray(:,2),'exp1');
  density=f(altitude);
  
  %plot(densityarray(:,1),densityarray(:,2))
  %hold on;
  %plot(f)

  r0=meanRadiusOfEarth+altitude; %% in m
  v=sqrt(mu/r0);
  MeanMotion=sqrt(mu/r0^3);          %% mean motion 
end

%{ 
https://ccmc.gsfc.nasa.gov/cgi-bin/modelweb/models/vitmo_model.cgi
VITMO ModelWeb Browser Results
NRLMSISE-00 model listing
Input parameters
year= 2000, month= 1, day= 1, hour= 1.50,
Time_type = Universal
Coordinate_type = Geographic
latitude= 0.00, longitude= 0.00, height= 100.00
Prof. parameters: start= 300.00 stop= 600.00 step= 50.00
Optional parametes: F10.7(daily) =not specified; F10.7(3-month avg) =not specified; ap(daily) = not specified
   Selected parameters are:
1 Height, km
2 Mass_density, g/cm-3
      1         2
  300.0 2.458E-14
  350.0 9.025E-15
  400.0 3.560E-15
  450.0 1.472E-15
  500.0 6.304E-16
  550.0 2.784E-16
  600.0 1.270E-16
%}