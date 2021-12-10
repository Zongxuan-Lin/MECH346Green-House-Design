close all;clear all;clc
%Cabin shaped greenhouse heat transfer numerical model
%%Structural Analysis,make sure it won't collpase in event of extreme
%%weather
%Brossard DataSnow load 3kpa, wind load 4kpa
%normal importance 
%10mm polycarbonate sheet density 0.35lb/sqft 1.70885kg/m^2
wat=1.70885*1*3*2^0.5;
was=1.70885*1*3*2;
woakroof=712*0.0254*0.1524*(4.5+3*2^0.5);
wroof=wat+woakroof
sd=wroof/(0.0254*0.1524*2)
roofloadtotal=sd+(3*3/2+4*3/2)*1000%much smaller than oak strengh the structure is rigid .
%uderdefined variables
W=3;    %Width of cabin
L=6;    %Length of cabin
W2=sqrt(2)*W/2;     %Roof WIdth
H1=1.5;     %wall height
H2=W/2;     %Roof Height
Qr=300;     %Radiation per hour, a total of 6 hours
Vwind=16;   %mean wind speed
tdb=-20;     %outside temperature
tsky=tdb-10;     %sky temperature
T2=17;      %Inside max temperature
T1=15;      %Inside minimum temperature
tAVG=(T1+T2)/2;     %average inside temperature 
tbt=-10;    %graound temperature below the platform
lb=20;       %number of layers of polyethyane films
la=1;   %Number of layers of polycarbonate board
%%
%Air properties
%Air at Edge@T1
rouin=(1.204-1.225)*(T1-15)/5+1.225; %Density
cpin=1007; %Specific Heat
Kin=(0.02514-0.02476)*(T1-15)/5+0.2476;
alpin=(2.074*10^(-5)-2.009*10^(-5))*(T1-15)/5+2.009*10^(-5);
miuin=(1.825*10^(-5)-1.802*10^(-5))*(T1-15)/5+1.802*10^(-5);
vin=(1.516*10^(-5)-1.47*10^(-5))*(T1-15)/5+1.47*10^(-5);

%Exterior Convection Air at -20
rouex=1.394;
cpex=1005;
Kex=0.02211;
alpex=1.578*10^(-5);
miuex=1.63*10^(-5);
vex=1.169*10^(-5);
prex=0.7408;

%Exterior Convection Air at -10
roubt=1.341;
cpbt=1006;
Kbt=0.02288;
alpbt=1.696*10^(-5);
miubt=1.68*10^(-5);
vbt=1.252*10^(-5);
prbt=0.7387;

% %Exterior Convection Air at 0
% roubt=1.292;
% cpbt=1006;
% Kbt=0.02364;
% alpbt=1.818*10^(-5);
% miubt=1.729*10^(-5);
% vbt=1.338*10^(-5);
% prbt=0.7362;
%%
%Construction material
%Polycarbonate Wall
ww=0.0103;  %panel thickness

K=0.0668;   %R=thickness/(K*A)%From a lab measuring the conductivity
n=1.59;      %index of reflection
em=0.95;     %emmisivity 
taostar=(1-((1-n)/(1+n))^2)^2/(1-((1-n)/(1+n))^4)^la;
eab=Qr*6*taostar;
%Conductive Resistance
RCconPC=ww/(K*L*H1)*la;
RBconPC=ww/(K*L*W2)*la;
RAconPC=ww/(K*(W*H1+W*H2/2))*la;

%polyethylene  Film 6mil,avag distance between film 1mm
ww1=0.0001524;
ww3=0.001;
K1=0.5;
n1=1.54;
em1=0.1;
taostar1=(1-((1-n1)/(1+n1))^2)^2/(1-((1-n1)/(1+n1))^4)^lb;
eab=eab*taostar1;
%Conductive Resistance
RCconPE=ww1/(K1*L*H1)*lb;
RBconPE=ww1/(K1*L*W2)*lb;
RAconPE=ww1/(K1*(W*H1+W*H2/2))*lb;
%conductive resistance of air in between films,using air property at
%-10,give a little overestimation
RCconCA=ww3/(Kbt*L*H1)*lb;
RBconCA=ww3/(Kbt*L*W2)*lb;
RAconCA=ww3/(Kbt*(W*H1+W*H2/2))*lb;

%Floor Conductive Resistance
%Assuming -10 of ari below the deck assume total conduction
tbelow=0.127; 
%oak wood board hardwood structure cp=1260 dens=721,
Kstr=0.17;
tstr=0.0256;
%%overall conductive thermal resistance
Rbl=tstr/(Kstr*W*L)+tbelow/(Kbt*W*L);
RAcond=RAconCA+RAconPE+RAconPC;
RBcond=RBconCA+RBconPE+RBconPC;
RCcond=RCconCA+RCconPE+RBconPC;
%%
%Neglect Convective heat transfer inside
if lb==0
    emasbl=0.95;
else
    emasbl=0.1;
end
surftemp=zeros(6,1);%case1 face a b c;case2 face a b c

%Case1 face a
t1a=tdb+10;t1b=tdb;t1c=(t1a+t1b)/2;Error=0.0001;A=W*H1+W*H2/2;
ha1a=convca1(t1a,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
ha1b=convca1(t1b,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
ha1c=convca1(t1c,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
fa=f(t1a,emasbl,tdb,ha1a,RAcond,tAVG,A);
fb=f(t1b,emasbl,tdb,ha1b,RAcond,tAVG,A);
fc=f(t1c,emasbl,tdb,ha1c,RAcond,tAVG,A);
while abs(fc)>Error
    if fc*fa<0
        t1b=t1c;
    else
        t1a=t1c;
    end
    t1c=(t1a+t1b)/2;
    ha1a=convca1(t1a,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
    ha1b=convca1(t1b,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
    ha1c=convca1(t1c,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
    fa=f(t1a,emasbl,tdb,ha1a,RAcond,tAVG,A);
    fb=f(t1b,emasbl,tdb,ha1b,RAcond,tAVG,A);
    fc=f(t1c,emasbl,tdb,ha1c,RAcond,tAVG,A);
end
surftemp(1)=t1c;

%Case1 face b
t1a=tdb+10;t1b=tdb;t1c=(t1a+t1b)/2;Error=0.0001;A=L*H1;
ha1a=convcb1(t1a,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,W2,L);
ha1b=convcb1(t1b,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,W2,L);
ha1c=convcb1(t1c,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,W2,L);
fa=f(t1a,emasbl,tdb,ha1a,RBcond,tAVG,A);
fb=f(t1b,emasbl,tdb,ha1b,RBcond,tAVG,A);
fc=f(t1c,emasbl,tdb,ha1c,RBcond,tAVG,A);
while abs(fc)>Error
    if fc*fa<0
        t1b=t1c;
    else
        t1a=t1c;
    end
    t1c=(t1a+t1b)/2;
    ha1a=convcb1(t1a,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,W2,L);
    ha1b=convcb1(t1b,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,W2,L);
    ha1c=convcb1(t1c,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,W2,L);
    fa=f(t1a,emasbl,tdb,ha1a,RBcond,tAVG,A);
    fb=f(t1b,emasbl,tdb,ha1b,RBcond,tAVG,A);
    fc=f(t1c,emasbl,tdb,ha1c,RBcond,tAVG,A);
end
surftemp(2)=t1c;

%Case1 face c
t1a=tdb+10;t1b=tdb;t1c=(t1a+t1b)/2;Error=0.0001;A=L*W2;
ha1a=convcc1(t1a,tdb,Kex,Kbt,tbt,vex,vbt,L,prex,prbt,Vwind,H1);
ha1b=convcc1(t1b,tdb,Kex,Kbt,tbt,vex,vbt,L,prex,prbt,Vwind,H1);
ha1c=convcc1(t1c,tdb,Kex,Kbt,tbt,vex,vbt,L,prex,prbt,Vwind,H1);
fa=f(t1a,emasbl,tdb,ha1a,RCcond,tAVG,A);
fb=f(t1b,emasbl,tdb,ha1b,RCcond,tAVG,A);
fc=f(t1c,emasbl,tdb,ha1c,RCcond,tAVG,A);
while abs(fc)>Error
    if fc*fa<0
        t1b=t1c;
    else
        t1a=t1c;
    end
    t1c=(t1a+t1b)/2;
    ha1a=convcc1(t1a,tdb,Kex,Kbt,tbt,vex,vbt,L,prex,prbt,Vwind,H1);
    ha1b=convcc1(t1b,tdb,Kex,Kbt,tbt,vex,vbt,L,prex,prbt,Vwind,H1);
    ha1c=convcc1(t1c,tdb,Kex,Kbt,tbt,vex,vbt,L,prex,prbt,Vwind,H1);
    fa=f(t1a,emasbl,tdb,ha1a,RCcond,tAVG,A);
    fb=f(t1b,emasbl,tdb,ha1b,RCcond,tAVG,A);
    fc=f(t1c,emasbl,tdb,ha1c,RCcond,tAVG,A);
end
surftemp(3)=t1c;

%Case2 face a
t1a=tdb+10;t1b=tdb;t1c=(t1a+t1b)/2;Error=0.0001;A=W*H1+W*H2/2;
ha1a=convca2(t1a,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
ha1b=convca2(t1b,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
ha1c=convca2(t1c,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
fa=f(t1a,emasbl,tdb,ha1a,RAcond,tAVG,A);
fb=f(t1b,emasbl,tdb,ha1b,RAcond,tAVG,A);
fc=f(t1c,emasbl,tdb,ha1c,RAcond,tAVG,A);
while abs(fc)>Error
    if fc*fa<0
        t1b=t1c;
    else
        t1a=t1c;
    end
    t1c=(t1a+t1b)/2;
    ha1a=convca2(t1a,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
    ha1b=convca2(t1b,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
    ha1c=convca2(t1c,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1);
    fa=f(t1a,emasbl,tdb,ha1a,RAcond,tAVG,A);
    fb=f(t1b,emasbl,tdb,ha1b,RAcond,tAVG,A);
    fc=f(t1c,emasbl,tdb,ha1c,RAcond,tAVG,A);
end
surftemp(4)=t1c;

%Case2 face b
t1a=tdb+10;t1b=tdb;t1c=(t1a+t1b)/2;Error=0.0001;A=W*H1+W*H2/2;
ha1a=convcb2(t1a,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,W2,L);
ha1b=convcb2(t1b,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,W2,L);
ha1c=convcb2(t1c,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,W2,L);
fa=f(t1a,emasbl,tdb,ha1a,RAcond,tAVG,A);
fb=f(t1b,emasbl,tdb,ha1b,RAcond,tAVG,A);
fc=f(t1c,emasbl,tdb,ha1c,RAcond,tAVG,A);
while abs(fc)>Error
    if fc*fa<0
        t1b=t1c;
    else
        t1a=t1c;
    end
    t1c=(t1a+t1b)/2;
    ha1a=convcb2(t1a,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,W2,L);
    ha1b=convcb2(t1b,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,W2,L);
    ha1c=convcb2(t1c,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,W2,L);
    fa=f(t1a,emasbl,tdb,ha1a,RAcond,tAVG,A);
    fb=f(t1b,emasbl,tdb,ha1b,RAcond,tAVG,A);
    fc=f(t1c,emasbl,tdb,ha1c,RAcond,tAVG,A);
end
surftemp(5)=t1c;

%Case2 face c
t1a=tdb+10;t1b=tdb;t1c=(t1a+t1b)/2;Error=0.0001;A=W*H1+W*H2/2;
ha1a=convcc2(t1a,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,H1,L);
ha1b=convcc2(t1b,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,H1,L);
ha1c=convcc2(t1c,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,H1,L);
fa=f(t1a,emasbl,tdb,ha1a,RAcond,tAVG,A);
fb=f(t1b,emasbl,tdb,ha1b,RAcond,tAVG,A);
fc=f(t1c,emasbl,tdb,ha1c,RAcond,tAVG,A);
while abs(fc)>Error
    if fc*fa<0
        t1b=t1c;
    else
        t1a=t1c;
    end
    t1c=(t1a+t1b)/2;
    ha1a=convcc2(t1a,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,H1,L);
    ha1b=convcc2(t1b,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,H1,L);
    ha1c=convcc2(t1c,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,H1,L);
    fa=f(t1a,emasbl,tdb,ha1a,RAcond,tAVG,A);
    fb=f(t1b,emasbl,tdb,ha1b,RAcond,tAVG,A);
    fc=f(t1c,emasbl,tdb,ha1c,RAcond,tAVG,A);
end
surftemp(6)=t1c;
surftemp2=zeros(6,1);%PC Wall thermal Radiation 
surftemp
surftemp2(1)=(surftemp(1)-tAVG)/RAcond*RAconCA+tAVG;
surftemp2(2)=(surftemp(2)-tAVG)/RBcond*RBconCA+tAVG;
surftemp2(3)=(surftemp(3)-tAVG)/RCcond*RCconCA+tAVG;
surftemp2(4)=(surftemp(4)-tAVG)/RAcond*RAconCA+tAVG;
surftemp2(5)=(surftemp(5)-tAVG)/RBcond*RBconCA+tAVG;
surftemp2(6)=(surftemp(6)-tAVG)/RCcond*RCconCA+tAVG;
Radinter=zeros(6,1);
Radinter(1)=(H1*W+W*H2/2)*5.67*10^(-8)*((surftemp2(1)+273.15)^4-(surftemp(1)+273.15)^4)/(1/em+(1/em1)*(2*lb-1)-lb);
Radinter(2)=(L*W2)*5.67*10^(-8)*((surftemp2(2)+273.15)^4-(surftemp(2)+273.15)^4)/(1/em+(1/em1)*(2*lb-1)-lb);
Radinter(3)=(L*H1)*5.67*10^(-8)*((surftemp2(3)+273.15)^4-(surftemp(3)+273.15)^4)/(1/em+(1/em1)*(2*lb-1)-lb);
Radinter(4)=(H1*W+W*H2/2)*5.67*10^(-8)*((surftemp2(4)+273.15)^4-(surftemp(4)+273.15)^4)/(1/em+(1/em1)*(2*lb-1)-lb);
Radinter(5)=(L*W2)*5.67*10^(-8)*((surftemp2(5)+273.15)^4-(surftemp(5)+273.15)^4)/(1/em+(1/em1)*(2*lb-1)-lb);
Radinter(6)=(L*H1)*5.67*10^(-8)*((surftemp2(6)+273.15)^4-(surftemp(6)+273.15)^4)/(1/em+(1/em1)*(2*lb-1)-lb);
Radinter%Too small and thus negligible
%%
%case1 face a b c;case2 face a b c
%Rbl=tstr/(Kstr*W*L)+tbelow/(Kbt*W*L);
% RAcond=RAconCA+RAconPE+RAconPC;
% RBcond=RBconCA+RBconPE+RBconPC;
% RCcond=RCconCA+RCconPE+RBconPC;
%case1 heat transfer rate day
powerin=eab/6*W*L
%case1 heat transfer rate 
Q1=(tAVG-tbt)/Rbl+(tAVG-surftemp(1))/RAcond*2+(tAVG-surftemp(2))/RBcond*2+(tAVG-surftemp(3))/RCcond*2+Radinter(1)+Radinter(2)+Radinter(3)
%case2 heat transfer rate 
Q2=(tAVG-tbt)/Rbl+(tAVG-surftemp(4))/RAcond*2+(tAVG-surftemp(5))/RBcond*2+(tAVG-surftemp(6))/RCcond*2+Radinter(4)+Radinter(5)+Radinter(6)
%Estimation of number of batteries and solar panel needed
Qavgtot=(Q1+Q2)/2*18/0.9/0.2;
Asp=Qavgtot/Qr/6
%%
%External Convection,Use Tdb as approximation since the outside surface
%temperature is quite close to the Tdb
%Case 1 wind parallel to door
%Six functions belwo calculate the inverse of convective thermal resistance
%h*A
function ha1=convca1(t,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1)
t1=(t+tdb)/2;
Ksf=(Kex-Kbt)*(t1-tbt)/(tdb-tbt)+Kbt;
Vsf=(vex-vbt)*(t1-tbt)/(tdb-tbt)+vbt;
Prsf=(prex-prbt)*(t1-tbt)/(tdb-tbt)+prbt;
Rea=Vwind*W/Vsf;
NumMesh=10;
ha1=0;
Reatop=zeros(NumMesh,1);%Mesh the roof triangular panel
Nutop=zeros(NumMesh,1);
for i=1: NumMesh
    Reatop(i)=Rea*((10-i)/10+0.05);
end
for i=1: NumMesh
    if Reatop(i)<5*10^5
        Nutop(i)=0.664*(Reatop(i))^0.5*Prsf^(1/3);
    else
        Nutop(i)=0.037*(Reatop(i))^0.8*Prsf^(1/3);
    end
    ha1=ha1+Nutop(i)*Ksf/(W*(10-i)/10+0.05)*(W*(10-i)/10+0.05)*H2/10;
end
if Rea<5*10^5
	Nua=0.664*(Rea)^0.5*Prsf^(1/3);
else
	Nua=0.037*(Rea)^0.8*Prsf^(1/3);
end
ha1=ha1+Nua*Ksf/(W)*W*H1;
end
function hb1=convcb1(t,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,W2,L)
t1=(t+tdb)/2;
Ksf=(Kex-Kbt)*(t1-tbt)/(tdb-tbt)+Kbt;
Vsf=(vex-vbt)*(t1-tbt)/(tdb-tbt)+vbt;
Prsf=(prex-prbt)*(t1-tbt)/(tdb-tbt)+prbt;
Reb=Vwind*W/Vsf;
Nub=0.258*(Reb)^0.588*Prsf^(1/3);
hb1=Nub*Ksf/(W)*W2*L;
end
function hc1=convcc1(t,tdb,Kex,Kbt,tbt,vex,vbt,L,prex,prbt,Vwind,H1)
t1=(t+tdb)/2;
Ksf=(Kex-Kbt)*(t1-tbt)/(tdb-tbt)+Kbt;
Vsf=(vex-vbt)*(t1-tbt)/(tdb-tbt)+vbt;
Prsf=(prex-prbt)*(t1-tbt)/(tdb-tbt)+prbt;
Rec=Vwind*L/Vsf;
Nuc=0.094*(Rec)^0.675*Prsf^(1/3);
hc1=Nuc*Ksf/(L)*L*H1;
end
%Case 2 wind in front
function ha2=convca2(t,tdb,Kex,Kbt,tbt,vex,vbt,W,prex,prbt,Vwind,H2,H1)
t1=(t+tdb)/2;
Ksf=(Kex-Kbt)*(t1-tbt)/(tdb-tbt)+Kbt;
Vsf=(vex-vbt)*(t1-tbt)/(tdb-tbt)+vbt;
Prsf=(prex-prbt)*(t1-tbt)/(tdb-tbt)+prbt;
Rea=Vwind*W/Vsf;
NumMesh=10;
ha2=0;
Reatop=zeros(NumMesh,1);%Mesh the roof triangular panel
Nutop=zeros(NumMesh,1);
for i=1: NumMesh
    Reatop(i)=Rea*((10-i)/10+0.05);
end
for i=1: NumMesh
    Nutop(i)=0.094*(Reatop(i))^0.675*Prsf^(1/3);
    ha2=ha2+Nutop(i)*Ksf/(W*(10-i)/10+0.05)*(W*(10-i)/10+0.05)*H2/10;
end
Nua=0.094*(Rea)^0.675*Prsf^(1/3);
ha2=ha2+Nua*Ksf/(W)*W*H1;
end
function hb2=convcb2(t,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,W2,L)
t1=(t+tdb)/2;
Ksf=(Kex-Kbt)*(t1-tbt)/(tdb-tbt)+Kbt;
Vsf=(vex-vbt)*(t1-tbt)/(tdb-tbt)+vbt;
Prsf=(prex-prbt)*(t1-tbt)/(tdb-tbt)+prbt;
Reb=Vwind*L/Vsf;
hb2=0;
if Reb<5*10^5
	Nub=0.664*(Reb)^0.5*Prsf^(1/3);
else
	Nub=0.037*(Reb)^0.8*Prsf^(1/3);
end
hb2=hb2+Nub*Ksf/(L)*W2*L;
end
function hc2=convcc2(t,tdb,Kex,Kbt,tbt,vex,vbt,prex,prbt,Vwind,H1,L)
t1=(t+tdb)/2;
Ksf=(Kex-Kbt)*(t1-tbt)/(tdb-tbt)+Kbt;
Vsf=(vex-vbt)*(t1-tbt)/(tdb-tbt)+vbt;
Prsf=(prex-prbt)*(t1-tbt)/(tdb-tbt)+prbt;
Rec=Vwind*L/Vsf;
hc2=0;
if Rec<5*10^5
	Nuc=0.664*(Rec)^0.5*Prsf^(1/3);
else
	Nuc=0.037*(Rec)^0.8*Prsf^(1/3);
end
hc2=hc2+Nuc*Ksf/(L)*H1*L;
end

%Function to find surface temperature 
function func = f(T,emasbl,tdb,hA,Rcond,tAVG,A)
    func = hA*(T-tdb)+5.67*10^(-8)*A*emasbl*((T+273.15)^4-(tdb+263.15)^4)-(tAVG-T)/(Rcond);
end

