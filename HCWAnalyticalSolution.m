%% analytical solution HCW
clc;clear all;close all;

altitude=340000;        %% in m
density=1e-8;
radiusOfEarth=6371000;  %% in m;

mu=3.986004418E14;      %% in m3?s?2
r0=radiusOfEarth+altitude; %% in m
omega=sqrt(mu/r0^3);

%% initial condition
x0=1;
y0=0;
z0=0;

u0=0;
v0=0;
w0=omega/2;

C1=u0/omega+2*z0
C2=w0/omega
C3=-3*z0-2*u0/omega
C4=x0-2*w0/omega
C5=v0/omega
C6=y0

t=0:60:2*pi/omega-pi/4/omega;

x=-3*C1*omega*t + 2*C2*cos(omega*t) - 2*C3*sin(omega*t) + C4;
y= C5*sin(omega*t) + C6*cos(omega*t);
z= 2*C1 + C2*sin(omega*t) + C3*cos(omega*t);

plot3(x,y,z)
hold on;
grid on;
vectarrow([1 0 0]);
text(1, 0, 0,"x",'HorizontalAlignment','left','FontSize',12);
hold on;
vectarrow([0 1 0]);
text(0, 1, 0,"y",'HorizontalAlignment','left','FontSize',12);
hold on;
vectarrow([0 0 1]);
text(0, 0, 1,"z",'HorizontalAlignment','left','FontSize',12);
axis equal;hold off;xlabel('x');ylabel('y');zlabel('z');axis([-1.55 1.55 -1.55 1.55 -1.55 1.55]);
text(x(end), y(end), z(end),"deputy",'HorizontalAlignment','left','FontSize',12);
text(0, 0, 0,"chief",'HorizontalAlignment','left','FontSize',12);
hold off;

function vectarrow(p1)
%Arrowline 3-D vector plot.
%   vectarrow(p0,p1) plots a line vector with arrow pointing from point p0
%   to point p1. The function can plot both 2D and 3D vector with arrow
%   depending on the dimension of the input
%
%   Example:
%       3D vector
%       p0 = [1 2 3];   % Coordinate of the first point p0
%       p1 = [4 5 6];   % Coordinate of the second point p1
%       vectarrow(p0,p1)
%
%       2D vector
%       p0 = [1 2];     % Coordinate of the first point p0
%       p1 = [4 5];     % Coordinate of the second point p1
%       vectarrow(p0,p1)
%
%   See also Vectline
%   Rentian Xiong 4-18-05
%   $Revision: 1.0
  
%% added by Jan
  p0 = [0 0 0];
%% end
  if max(size(p0))==3
      if max(size(p1))==3
          x0 = p0(1);
          y0 = p0(2);
          z0 = p0(3);
          x1 = p1(1);
          y1 = p1(2);
          z1 = p1(3);
          plot3([x0;x1],[y0;y1],[z0;z1]);   % Draw a line between p0 and p1
          
          p = p1-p0;
          alpha = 0.1;  % Size of arrow head relative to the length of the vector
          beta = 0.1;  % Width of the base of the arrow head relative to the length
          
          hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
          hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
          hw = [z1-alpha*p(3);z1;z1-alpha*p(3)];
          
          hold on
          plot3(hu(:),hv(:),hw(:))  % Plot arrow head
          grid on
          xlabel('x')
          ylabel('y')
          zlabel('z')
          hold off
      else
          error('p0 and p1 must have the same dimension')
      end
  elseif max(size(p0))==2
      if max(size(p1))==2
          x0 = p0(1);
          y0 = p0(2);
          x1 = p1(1);
          y1 = p1(2);
          plot([x0;x1],[y0;y1]);   % Draw a line between p0 and p1
          
          p = p1-p0;
          alpha = 0.1;  % Size of arrow head relative to the length of the vector
          beta = 0.1;  % Width of the base of the arrow head relative to the length
          
          hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
          hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
          
          hold on
          plot(hu(:),hv(:))  % Plot arrow head
          grid on
          xlabel('x')
          ylabel('y')
          hold off
      else
          error('p0 and p1 must have the same dimension')
      end
  else
      error('this function only accepts 2D or 3D vector')
  end
end