clear all;clc;close all;
%% location of the TLE files, they shall be zipped text files in twolines, i.e. with a name line
%% if you use the CelesTrak TLE retriever, then this setting will work
tleFiles=dir(strcat(getenv('USERPROFILE'),'\Documents\My TLEs\2019\Full Catalog\*'));
%% which satellites to be processed, use NORAD ID
catalogueID=[43765,43794,43799]; %% hawk-a,hawk-b,hawk-c
%catalogueID=[43196,43197]; %% GOMX-4B GOMX-4A


%% initialize arrays
epochTime=zeros( size(catalogueID,2) ,1);
RAAN=zeros( size(catalogueID,2) ,1);
a=zeros( size(catalogueID,2) ,1);
w=zeros( size(catalogueID,2) ,1);
dn=zeros( size(catalogueID,2) ,1);
Incl=zeros( size(catalogueID,2) ,1);

epochTime2=zeros( size(catalogueID,2) ,1);
RAAN2=zeros( size(catalogueID,2) ,1);
a2=zeros( size(catalogueID,2) ,1);
w2=zeros( size(catalogueID,2) ,1);
dn2=zeros( size(catalogueID,2) ,1);
Incl2=zeros( size(catalogueID,2) ,1);


for i=3:size(tleFiles,1) %% %% cycle overall all files, first two entries of tleFiles are '.' and '..'. They shall be skipped
  %% copy file here
  copyfile(strcat(tleFiles(i).folder,'\',tleFiles(i).name), 'temp/');
  %% unzip file
  unzippedfilename=unzip(strcat('temp/',tleFiles(i).name),'temp/');
  %% readfile into variable
  [epochTimeN,RAANN,aN,wN,InclN,foundno]=readtle( unzippedfilename{1},catalogueID);
  %% add TLE of catalogueID in arrary
  if size(epochTimeN,2)==1 && size(epochTimeN,1)==size(catalogueID,2) && sum(foundno)~=0
    dnN=zeros( size(catalogueID,2) ,1);
    for j=1:size(epochTimeN,1)
        if epochTimeN(j)~=0
          datestringN=strcat('20',num2str(epochTimeN(j)));
          dnN(j)=datenum(str2num(datestringN(1:4)),1,str2num(datestringN(5:10)));
        else
          dnN(j)=0;
        end
    end
    %% add new entries
    epochTime=[epochTime epochTimeN];
    dn=[dn dnN];
    RAAN=[RAAN RAANN];
    a=[a aN];
    w=[w wN];
    Incl=[Incl InclN];
    %for j=1:size(catalogueID,2)
    %  sat2(j).epochTime=[sat2(j).epochTime epochTimeN(j)];
    %end
  elseif size(epochTimeN,2)>1 && size(epochTimeN,1)==size(catalogueID,2) && sum(foundno)~=0
    dnN2=zeros(size(epochTimeN));   
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
    %% add new entries
    epochTime2=[epochTime2 epochTimeN];
    dn2=[dn2 dnN2];
    RAAN2=[RAAN2 RAANN];
    a2=[a2 aN];
    w2=[w2 wN];
    Incl2=[Incl2 InclN];
  elseif size(epochTimeN,2)==1 && size(epochTimeN,1)==size(catalogueID,2) && sum(foundno)~=1
    ;
  else
    fprintf('\n error');
    input('');
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
Incl(:,1)=[];

epochTime2(:,1)=[];
RAAN2(:,1)=[];
a2(:,1)=[];
w2(:,1)=[];
dn2(:,1)=[];
Incl2(:,1)=[];


for i=1:size(catalogueID,2)
  sat(i).epochTime=epochTime(i,:);
  sat(i).RAAN=RAAN(i,:);
  sat(i).a=a(i,:);
  sat(i).w=w(i,:);
  sat(i).dn=dn(i,:);
  sat(i).Incl=Incl(i,:);
  
  sat(i).epochTime=[sat(i).epochTime epochTime2(i,epochTime2(i,:)~=0)];
  sat(i).RAAN=[sat(i).RAAN RAAN2(i,RAAN2(i,:)~=0)];
  sat(i).a=[sat(i).a a2(i,a2(i,:)~=0)];
  sat(i).w=[sat(i).w w2(i,w2(i,:)~=0)];
  sat(i).dn=[sat(i).dn dn2(i,dn2(i,:)~=0)];
  sat(i).Incl=[sat(i).Incl Incl2(i,Incl2(i,:)~=0)];

  %% define time and order arrary accordingly
  [sat(i).epochTime idx]=sort(sat(i).epochTime);
  sat(i).RAAN=sat(i).RAAN(idx);
  sat(i).a=sat(i).a(idx);
  sat(i).w=sat(i).w(idx);
  sat(i).dn=sat(i).dn(idx);
  sat(i).Incl=sat(i).Incl(idx);
  
  %% interpolate on first satellite' time instances
  if i~=1
    [TEMP, idx]=unique(sat(i).dn);
    sat(i).RAANInt     =interp1(TEMP,sat(i).RAAN(idx),sat(1).dn);
    sat(i).aInt        =interp1(TEMP,sat(i).a(idx),sat(1).dn);
    sat(i).wInt        =interp1(TEMP,sat(i).w(idx),sat(1).dn);
    sat(i).InclInt     =interp1(TEMP,sat(i).Incl(idx),sat(1).dn);
  end  
end


%% absolut plot argument of perigee, RAAN, excentrity
figure
  set(gcf, 'Position',  [50, 50, 1500, 500]);
  subplot(2,4,1)
    for i=1:size(catalogueID,2)
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).RAAN);hold on;
    end
    ylabel('RAAN [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    legend;
  subplot(2,4,2)
    for i=1:size(catalogueID,2)
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).Incl);hold on;
    end
    ylabel('incl [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
  subplot(2,4,3)
    for i=1:size(catalogueID,2)    
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).a-6371000);hold on;
    end
    ylabel('alt [m]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    %axis([0 sat(1).dn(end)-sat(1).dn(1) mean(sat(1).a-6371000)-0.001*std(sat(1).a-6371000) mean(sat(1).a-6371000)+0.06*std(sat(1).a-6371000)])
  subplot(2,4,4)
    for i=1:size(catalogueID,2)  
      plot(sat(i).dn(:)-sat(i).dn(1),sat(i).w);hold on;
    end
    ylabel('a o peri [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );

%% relative plot argument of perigee, RAAN, excentrity (relative to first satellite)
  subplot(2,4,5)
    for i=2:size(catalogueID,2)
      plot(sat(1).dn-sat(1).dn(1),sat(i).RAANInt-sat(1).RAAN);hold on;
    end
    ylabel('rel RAAN [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.001*std(sat(1).RAAN) 0.001*std(sat(1).RAAN)]);
    grid on; legend;   
  subplot(2,4,6)
    for i=2:size(catalogueID,2)
      plot(sat(1).dn-sat(1).dn(1),sat(i).InclInt-sat(1).Incl);hold on;
    end
    ylabel('rel incl [deg]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.001*std(sat(1).Incl) 0.001*std(sat(1).Incl)]);
    grid on;  
  subplot(2,4,7)
    for i=2:size(catalogueID,2)    
      plot(sat(1).dn-sat(1).dn(1),sat(i).aInt-sat(1).a);hold on;
    end
    ylabel('rel alt [m]'); xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    %axis([0 sat(1).dn(end)-sat(1).dn(1) -0.002*std(sat(1).a) 0.002*std(sat(1).a)]);
    grid on;
  subplot(2,4,8)
    for i=2:size(catalogueID,2)  
      plot(sat(1).dn-sat(1).dn(1),sat(i).wInt-sat(1).w);hold on;
    end
    ylabel('rel a o peri[deg]');xlabel(strcat('time from',{' '},datestr(datetime( sat(1).dn(1),'ConvertFrom','datenum') ),{' '},'[d]') );
    axis([0 sat(1).dn(end)-sat(1).dn(1) -0.01*std(sat(1).w) 0.01*std(sat(1).w)]);
    grid on; 

function [epochTime,RAAN,a,w,Incl,foundno]=readtle(file, catalog)

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
  foundno=zeros(size(catalog,2),1);
  
  epochTime=zeros(size(catalog,2),1);
  Incl=zeros(size(catalog,2),1);
  RAAN=zeros(size(catalog,2),1);
  ecc=zeros(size(catalog,2),1);
  w=zeros(size(catalog,2),1);
  M=zeros(size(catalog,2),1);
  n=zeros(size(catalog,2),1);
  T=zeros(size(catalog,2),1);
  a=zeros(size(catalog,2),1);
  b=zeros(size(catalog,2),1);
  
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
