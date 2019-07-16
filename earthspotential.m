clc;clear all;close all;


L=0/180*pi:pi/10:360/180*pi;
L=L;
mu=3.986004418E14;      %% in m3?s?2
radiusOfEarth=6371000;
J=[1 0.00108263 -0.00000254 -0.00000161];
A=zeros(size(L,2),1);
V=zeros(size(L,2),4);

 for n=2:4
    A=J(n)*(radiusOfEarth/radiusOfEarth)^n*legendreP(n,sin(L'));
    V(:,n)=mu/radiusOfEarth*(1-A);
 end
 
figure

subplot(3,1,1)
polarplot(L,V(:,2));hold on;rlim([6.24e7 6.266e7]);legend(strcat('J',int2str(2)));thetalim([-90 90]);
subplot(3,1,2)
polarplot(L,V(:,3));hold on;rlim([6.2564e7 6.2566e7]);legend(strcat('J',int2str(3)));thetalim([-90 90]);
subplot(3,1,3)
polarplot(L,V(:,4));hold on;rlim([6.2564e7 6.2566e7]);legend(strcat('J',int2str(4)));thetalim([-90 90]);
