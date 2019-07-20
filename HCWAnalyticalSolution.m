%% analytical solution HCW
clc;clear all;close all;

altitude=340000;        %% in m
density=1e-8;
radiusOfEarth=6371000;  %% in m;

mu=3.986004418E14;      %% in m3?s?2
r0=radiusOfEarth+altitude; %% in m
omega=sqrt(mu/r0^3);
%omega=2*pi;

%% initial condition
x0=1;
y0=0;
z0=0;

u0=0;
v0=0;
w0=omega/2;

C1=u0/omega+2*z0;
C2=w0/omega;
C3=-3*z0-2*u0/omega;
C4=x0-2*w0/omega;
C5=v0/omega;
C6=y0;
t=0:(2*pi/omega)/100:(2*pi/omega);%T;%-T/8;


x=-3*C1*omega*t + 2*C2*cos(omega*t) - 2*C3*sin(omega*t) + C4;
y= C5*sin(omega*t) + C6*cos(omega*t);
z= 2*C1 + C2*sin(omega*t) + C3*cos(omega*t);

figure
  subplot(3,1,1)
  plot(t,x)
  subplot(3,1,2)
  plot(t,y)
  subplot(3,1,3)
  plot(t,z)

figure
  plot3(x,y,z)
  hold on; grid on;
  vectarrow([0 0 0],[1 0 0]);
  text(1, 0, 0,"x",'HorizontalAlignment','left','FontSize',12);
  hold on;
  vectarrow([0 0 0],[0 1 0]);
  text(0, 1, 0,"y",'HorizontalAlignment','left','FontSize',12);
  hold on;
  vectarrow([0 0 0],[0 0 1]);
  text(0, 0, 1,"z",'HorizontalAlignment','left','FontSize',12);
  axis equal;hold off;xlabel('x');ylabel('y');zlabel('z');%axis([-1.55 1.55 -1.55 1.55 -1.55 1.55]);
  text(x(end), y(end), z(end),"deputy",'HorizontalAlignment','left','FontSize',12);
  text(0, 0, 0,"chief",'HorizontalAlignment','left','FontSize',12);
  hold off;
