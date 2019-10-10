clear all;clc;close all;
tleFiles=dir('C:\Users\jan.thoemel\Documents\My TLEs\2019\Full Catalog\*');
catalogueID=[43765,43794,43799]; %% hawk-a,hawk-b,hawk-c

epochTime=zeros( size(catalogueID,2) ,1);
Omega=zeros( size(catalogueID,2) ,1);
a=zeros( size(catalogueID,2) ,1);
w=zeros( size(catalogueID,2) ,1);

for i=3:size(tleFiles,1) %% first two entries of tleFiles are '.' and '..'. They shall be skipped
  %% copy file here
  copyfile(strcat(tleFiles(i).folder,'\',tleFiles(i).name), 'temp/');
  %% unzip file
  unzippedfilename=unzip(strcat('temp/',tleFiles(i).name),'temp/');
  %% readfile into variable
  [epochTimeN,OmegaN,aN,wN]=readtle( unzippedfilename{1},catalogueID);
  epochTime=[epochTime epochTimeN'];
  Omega=[Omega OmegaN'];
  a=[a aN'];
  w=[w wN'];
  %% delete zip and unzipped file
  delete(strcat('temp/',tleFiles(i).name));
  delete(unzippedfilename{1});
  %%add TLE of Hawkeye in arrary
  
end

%% remove first column
epochTime(:,1)=[];
Omega(:,1)=[];
a(:,1)=[];
w(:,1)=[];


%% define time and order arrary accordingly



%% plot argument of perigee, RAAN, excentrity

  subplot(3,1,1)
    for i=1:size(catalogueID,2)
      plot(epochTime(i,:),Omega(i,:));hold on;
    end
  subplot(3,1,2)
    for i=1:size(catalogueID,2)    
      plot(epochTime(i,:),a(i,:));hold on;
    end
  subplot(3,1,3)
    for i=1:size(catalogueID,2)  
      plot(epochTime(i,:),w(i,:));hold on;
    end


function [epochTime,Omega,a,w]=readtle(file, catalog)

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
  
  %A0 = fgetl(fd);
  A1 = fgetl(fd);
  A2 = fgetl(fd);
  while ischar(A2)
    n = n + 1;
    satnum = str2num(A1(3:7));
    if isempty(catalog) || ismember(satnum, catalog)
      satcount=find(satnum==catalog);
      %fprintf('%s\n', repmat('-',1,50));
      %fprintf('Satellite: %s\n', A0)
      assert(chksum(A1), 'Checksum failure on line 1')
      assert(chksum(A2), 'Checksum failure on line 2')
      %fprintf('Catalog Number: %d\n', satnum)
      %fprintf('Epoch time: %s\n', A1(19:32)) % YYDDD.DDDDDDDD
      A1(19:32);
      epochTime(satcount)=str2double(A1(19:32));
      Incl(satcount) = str2num(A2(9:16));
      %fprintf('Inclination: %f deg\n', Incl)
      Omega(satcount) = str2num(A2(18:25));
      %fprintf('RA of ascending node: %f deg\n', Omega)
      ecc(satcount) = str2num(['.' A2(27:33)]);
      %fprintf('Eccentricity: %f\n', ecc)
      w(satcount) = str2num(A2(35:42));
      %fprintf('Arg of perigee: %f deg\n', w)
      M(satcount) = str2num(A2(44:51));
      %fprintf('Mean anomaly: %f deg\n', M)
      n(satcount) = str2num(A2(53:63));
      %fprintf('Mean motion: %f rev/day\n', n)
      T(satcount) = 86400/n(satcount);
      %fprintf('Period of rev: %.0f s/rev\n', T)
      a(satcount) = ((T(satcount)/(2*pi))^2*398.6e12)^(1/3);
      %fprintf('Semi-major axis: %.0f meters\n', a)
      b(satcount) = a(satcount)*sqrt(1-ecc(satcount)^2);
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
