
%% for testing
%OE=readTLE(fname);
%whichconst=0;
%opsmode=0;
%satrec=x;
%epoch=OE(1);
%xbstar=0;
%xndot=0;
%xnddot=0;
%xecco=OE(3);
%xargpo=OE(6);
%xinclo=OE(4);
%xmo=OE(7);
%xno_kozai=OE(9);
%xnodeo=OE(5);
%     satn        - satellite number
%     bstar       - sgp4 type drag coefficient              kg/m2er
%     ecco        - eccentricity
%     epoch       - epoch time in days from jan 0, 1950. 0 hr
%     argpo       - argument of perigee (output if ds)
%     inclo       - inclination
%     mo          - mean anomaly (output if ds)
%     no          - mean motion
%     nodeo      - right ascension of ascending node

%1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753
%2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667     0.00      4320.0        360.00

%fname='testtle.tle'

%OE=readTLE(fname);
%fprintf('epoch [day]           a [km]            e           inc [deg]        RAAN [deg]            w[deg]          M [deg]        Rev  MM \n ')
%OE
    
%sgp4init2(whichconst, opsmode, satrec, epoch, xbstar, xndot, xnddot,xecco, xargpo, xinclo, xmo, xno_kozai, xnodeo)
%% end for testing
clear all;clc;
longstr1='1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753';
longstr2='2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667';
typerun='c';
typeinput='d';
opsmode='i';
whichconst=84;
%     epoch       - epoch time in days from jan 0, 1950. 0 hr
epoch=365*69+100;
[startmfe, stopmfe, deltamin, satrec] = twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst);
[satrec] = sgp4init(whichconst, opsmode, satrec, epoch, satrec.bstar, satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai, satrec.nodeo);
tsince=0;
while tsince<stopmfe-startmfe
    tsince=tsince+10;
    [satrec, r, v] = sgp4(satrec, tsince);
    fprintf('%f %f %f\n',r(1),r(2),r(3));
    [latgc,latgd,lon,hellp] = ijk2lle ( r, jd );
end
