clear all;clc;close all;
tleFiles=dir('C:\Users\jan.thoemel\Documents\My TLEs\2019\Full Catalog\*');
catalogueID=[43765,43794,43799]; %% hawk-a,hawk-b,hawk-c

epochTime=zeros( size(catalogueID,2) ,1);
RAAN=zeros( size(catalogueID,2) ,1);
a=zeros( size(catalogueID,2) ,1);
w=zeros( size(catalogueID,2) ,1);
dn=zeros( size(catalogueID,2) ,1);

epochTime2=zeros( size(catalogueID,2) ,1);
RAAN2=zeros( size(catalogueID,2) ,1);
a2=zeros( size(catalogueID,2) ,1);
w2=zeros( size(catalogueID,2) ,1);
dn2=zeros( size(catalogueID,2) ,1);



for i=3:size(tleFiles,1) %% first two entries of tleFiles are '.' and '..'. They shall be skipped
  %% copy file here
  copyfile(strcat(tleFiles(i).folder,'\',tleFiles(i).name), 'temp/');
  %% unzip file
  unzippedfilename=unzip(strcat('temp/',tleFiles(i).name),'temp/');
  %% readfile into variable
  [epochTimeN,RAANN,aN,wN]=readtle( unzippedfilename{1},catalogueID);
  
  %% add TLE of catalogueID in arrary
  if size(epochTimeN,2)==1 && size(epochTimeN,1)==size(catalogueID,2)
    dnN=zeros( size(catalogueID,2) ,1);
    for j=1:size(epochTimeN,1)
      datestringN=strcat('20',num2str(epochTimeN(j)));
      dnN(j)=datenum(str2num(datestringN(1:4)),1,str2num(datestringN(5:10)));
    end
    epochTime=[epochTime epochTimeN];
    dn=[dn dnN];
    RAAN=[RAAN RAANN];
    a=[a aN];
    w=[w wN];
  else %% do nothing
    dnN2=zeros(size(epochTimeN));   
    %size(RAANN)
    %size(epochTime)
    %size(epochTimeN)
    %size(dnN)

    for j=1:size(epochTimeN,2)
      for k=1:size(epochTimeN,1)
        if epochTimeN(k,j)~=0
          datestringN=strcat('20',num2str(epochTimeN(k,j)));
          dnN2(k,j)=datenum(str2num(datestringN(1:4)),1,str2num(datestringN(5:end)));
        else
          dnN2(k,j)=0;
        end
      end
    end    
    epochTime2=[epochTime2 epochTimeN];
    dn2=[dn2 dnN2];
    RAAN2=[RAAN2 RAANN];
    a2=[a2 aN];
    w2=[w2 wN];    
  end
  %% delete zip and unzipped file
  delete(strcat('temp/',tleFiles(i).name));
  delete(unzippedfilename{1});  
  
  if i>3
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
  end
  fprintf('number of loaded TLEs %4.0f/%4.0f (%4.0f%%)', i, size(tleFiles,1),double(i)/double(size(tleFiles,1))*100 )
  pause(0.001);
  
end

%% remove first column
epochTime(:,1)=[];
RAAN(:,1)=[];
a(:,1)=[];
w(:,1)=[];
dn(:,1)=[];

epochTime2(:,1)=[];
RAAN2(:,1)=[];
a2(:,1)=[];
w2(:,1)=[];
dn2(:,1)=[];


for i=1:size(catalogueID,2)
  sat(i).epochTime=epochTime(i,:);
  sat(i).RAAN=RAAN(i,:);
  sat(i).a=a(i,:);
  sat(i).w=w(i,:);
  sat(i).dn=dn(i,:);
  sat(i).epochTimeInt=9;
  
  sat(i).epochTime=[sat(i).epochTime epochTime2(i,epochTime2(i,:)~=0)];
  sat(i).RAAN=[sat(i).RAAN RAAN2(i,RAAN2(i,:)~=0)];
  sat(i).a=[sat(i).a a2(i,a2(i,:)~=0)];
  sat(i).w=[sat(i).w w2(i,w2(i,:)~=0)];
  sat(i).dn=[sat(i).dn dn2(i,dn2(i,:)~=0)];
  
  [sat(i).epochTime idx]=sort(sat(i).epochTime);
  sat(i).RAAN=sat(i).RAAN(idx);
  sat(i).a=sat(i).a(idx);
  sat(i).w=sat(i).w(idx);
  sat(i).dn=sat(i).dn(idx);
  
  if i~=1
    [TEMP, idx]=unique(sat(i).dn);
    sat(i).RAANInt     =interp1(TEMP,sat(i).RAAN(idx),sat(1).dn);
    sat(i).aInt        =interp1(TEMP,sat(i).a(idx),sat(1).dn);
    sat(i).wInt        =interp1(TEMP,sat(i).w(idx),sat(1).dn);
  end  
end

%epochTime3(1,:)=epochTime2(1,epochTime2(1,:)~=0);

%% define time and order arrary accordingly


%% absolut plot argument of perigee, RAAN, excentrity
  subplot(3,1,1)
    for i=1:size(catalogueID,2)
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).RAAN);hold on;
    end
    ylabel('RAAN [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
  subplot(3,1,2)
    for i=1:size(catalogueID,2)    
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).a-6371000);hold on;
    end
    ylabel('altitude [m]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) mean(sat(1).a-6371000)-0.001*std(sat(1).a-6371000) mean(sat(1).a-6371000)+0.06*std(sat(1).a-6371000)])
  subplot(3,1,3)
    for i=1:size(catalogueID,2)  
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).w);hold on;
    end
    ylabel('arg of perigee [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );

%% relative plot argument of perigee, RAAN, excentrity (relative to first satellite)
figure
  subplot(3,1,1)
    for i=2:size(catalogueID,2)
      plot(sat(1).dn-sat(1).dn(1),sat(i).RAANInt-sat(1).RAAN);hold on;
    end
    ylabel('rel. RAAN [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.001*std(sat(1).RAAN) 0.001*std(sat(1).RAAN)]);grid on; legend;   
    subplot(3,1,2)
    for i=2:size(catalogueID,2)    
      plot(sat(1).dn-sat(1).dn(1),sat(i).aInt-sat(1).a);hold on;
    end
    ylabel('rel. altitude [m]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.002*std(sat(1).a) 0.002*std(sat(1).a)]);grid on;    
  subplot(3,1,3)
    for i=2:size(catalogueID,2)  
      plot(sat(1).dn-sat(1).dn(1),sat(i).wInt-sat(1).w);hold on;
    end
    ylabel('rel. arg of perigee [deg]');xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.01*std(sat(1).w) 0.01*std(sat(1).w)]);grid on; 


function [epochTime,RAAN,a,w]=readtle(file, catalog)

% READTLE Read satellite ephemeris data from a NORAD two-line element (TLE) file.
%
% INPUTS:
%   file    - Path to any standard two-line element file.
%   catalog - Optional array of NORAD catalog numbers for the satellites of
%             interest. The default action is to display data from every
%             satellite in the file.
%
% Brett Pantalone
% North Carolina State University
% Department of Electrical and Computer Engineering
% Optical Sensing Laboratory
% mailto:bapantal@ncsu.edu
% http://research.ece.ncsu.edu/osl/
  if nargin < 2
    catalog = [];
  end
  fd = fopen(file,'r');
  if fd < 0, fd = fopen([file '.tle'],'r'); end
  assert(fd > 0,['Can''t open file ' file ' for reading.'])
  n = 0;
  foundno=[0 0 0];

  epochTime=zeros(3,1);
  Incl=zeros(3,1);
  RAAN=zeros(3,1);
  ecc=zeros(3,1);
  w=zeros(3,1);
  M=zeros(3,1);
  n=zeros(3,1);
  T=zeros(3,1);
  a=zeros(3,1);
  b=zeros(3,1);
  
  %A0 = fgetl(fd);
  A1 = fgetl(fd);
  A2 = fgetl(fd);
  while ischar(A2)
    n = n + 1;
    satnum = str2num(A1(3:7));
    if isempty(catalog) || ismember(satnum, catalog)
      satcount=find(satnum==catalog);
      foundno(satcount)=foundno(satcount)+1;
      %fprintf('%s\n', repmat('-',1,50));
      %fprintf('Satellite: %s\n', A0)
      assert(chksum(A1), 'Checksum failure on line 1')
      assert(chksum(A2), 'Checksum failure on line 2')
      %fprintf('Catalog Number: %d\n', satnum)
      %fprintf('Epoch time: %s\n', A1(19:32)) % YYDDD.DDDDDDDD
      A1(19:32);
      epochTime(satcount,foundno(satcount))=str2double(A1(19:32));
      Incl(satcount,foundno(satcount)) = str2num(A2(9:16));
      %fprintf('Inclination: %f deg\n', Incl)
      RAAN(satcount,foundno(satcount)) = str2num(A2(18:25));
      %fprintf('RA of ascending node: %f deg\n', RAAN)
      ecc(satcount,foundno(satcount)) = str2num(['.' A2(27:33)]);
      %fprintf('Eccentricity: %f\n', ecc)
      w(satcount,foundno(satcount)) = str2num(A2(35:42));
      %fprintf('Arg of perigee: %f deg\n', w)
      M(satcount,foundno(satcount)) = str2num(A2(44:51));
      %fprintf('Mean anomaly: %f deg\n', M)
      n(satcount,foundno(satcount)) = str2num(A2(53:63));
      %fprintf('Mean motion: %f rev/day\n', n)
      T(satcount,foundno(satcount)) = 86400/n(satcount,foundno(satcount));
      %fprintf('Period of rev: %.0f s/rev\n', T)
      a(satcount,foundno(satcount)) = ((T(satcount,foundno(satcount))/(2*pi))^2*398.6e12)^(1/3);
      %fprintf('Semi-major axis: %.0f meters\n', a)
      b(satcount,foundno(satcount)) = a(satcount,foundno(satcount))*sqrt(1-ecc(satcount,foundno(satcount))^2);
      %fprintf('Semi-minor axis: %.0f meters\n', b)
    end
    %A0 = fgetl(fd);
    A1 = fgetl(fd);
    A2 = fgetl(fd);
  end
  fclose(fd);
end
%%
% Checksum (Modulo 10)
% Letters, blanks, periods, plus signs = 0; minus signs = 1
function result = chksum(str)
  result = false; c = 0;
  
  for k = 1:68
    if str(k) > '0' && str(k) <= '9'
      c = c + str(k) - 48;
    elseif str(k) == '-'
      c = c + 1;
    end
  end
  if mod(c,10) == str(69) - 48
    result = true;
  end
  
end
