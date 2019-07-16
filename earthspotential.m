clc;clear all;close all;


L=0/180*pi:pi/10:360/180*pi;
L=L;
mu=3.986004418E14;      %% in m3?s?2
radiusOfEarth=6378000;%%equator
J=[0 0.00108263 -0.00000254 -0.00000161];
A=zeros(size(L,2),1);
V=zeros(size(L,2),4);

 for n=1:4
    A=J(n)*(radiusOfEarth/radiusOfEarth)^n*legendreP(n,sin(L'));
    V(:,n)=mu/radiusOfEarth*(1-A);
 end
 
figure

subplot(2,2,1)
polarplot(L,V(:,1));rticks([]);rlim([5.23e7 7.266e7]);legend(strcat('J',int2str(1)));thetalim([-90 90]);
subplot(2,2,2)
polarplot(L,V(:,2));rticks([]);rlim([6.23e7 6.266e7]);legend(strcat('J',int2str(2)));thetalim([-90 90]);
subplot(2,2,3)
polarplot(L,V(:,3));rticks([]);rlim([6.2495e7 6.2497e7]);legend(strcat('J',int2str(3)));thetalim([-90 90]);
subplot(2,2,4)
polarplot(L,V(:,4));rticks([]);rlim([6.2495e7 6.24971e7]);legend(strcat('J',int2str(4)));thetalim([-90 90]);
