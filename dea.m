clc;
clear all;
close all;
%%
%Taking inputs 
Power=input('Enter the Rated Power of the motor: ');
choice=input('1. HP  2.kW: ');
if(choice==1)
        Po=Power*0.746;
else
        Po=Power;
end
Vin=input('Enter the rated voltage: ');
N=input('Enter the rated speed of the motor(in rpm): ');
Ns = N/60;%in rps
f=input('Enter the rated frequency of the motor: ');
Bav=input('Enter the specific magnetic loading: ');
pf=input('Enter the power factor of the motor: ');
eff=input('Enter the target efficiency of the motor: ');
ac=input('Enter the specific electrical loading: ');
choice2=input('\n1. Delta Connected 2.Star Connected: ');
%%
%Main dimension calculation
phase=3;%assumed
Q = Po/(eff*pf);%rating in kVA
Kw = .9;%assumed 
fprintf('\n SQUIRREL CAGE IND MOTOR DESIGN DATA');
fprintf('\n ——————————————------------');
Co = 1.11*pi*pi*Bav*ac*Kw*(10^-3);
%number of poles
pole = (2*f)/Ns;
%Q = Co*Ns*D2L
fprintf('\nInput power or rating power = ');
disp(Q);
D2L = Q/(Co*Ns);
%for good overall design consider L/tow = 1
% therefore L*pole/pi*D=1 or L*pole=Pi*D
D3 = (Q*pole)/(Co*Ns*pi);
D = D3^(1/3);
fprintf('\n MAIN DIMENSIONS ');
fprintf('\n --------------- ');
fprintf('\nHence Diameter D = ');
disp(D);
fprintf('\nHence Length L = ');
L = pi*D/(pole);
disp(L);
%%
%Other Characteristics Calculation
Va = pi*D*Ns;%peripheral speed Va
fprintf('\n CERTAIN OTHER CHARS OF THE MOTOR: ');
fprintf('\n -------------------------------- ');
fprintf('\nPeripheral speed = ');
disp(Va);
if(Va<30)
fprintf('\nAs Peripheral speed is less than 30m/secs so dimensions are permissable');
else
fprintf('\nAs Peripheral speed is not less than 30m/sec the dimensions are not mpermissabel. But still the dimensions will be');
end
phim = Bav*pi*D*L/pole;
fprintf('\nFlux density phim = ');
disp(phim);
fprintf('\n--------------------------------------------------------------');
%%
%Stator Characterisitics calculation
fprintf('\n STATOR CHARACTERISTICS ');
fprintf('\n ---------------------- ');
Ts = Vin/(4.44*f*phim*Kw);%number of stator turns Ts
Ts=ceil(Ts);
%------------------NUMBER OF STATOR SLOTS SELECTION----------------------
p=pole;
q=[2:1:20];
Ss=3*p*q;
yss=(pi*D)./Ss;
reqd=0;
for i=1:20
    if(yss(i)>0.015&yss(i)<0.025)
        reqd=i;
        break;
    end
end
yss=yss(reqd);
Ss=round(Ss(reqd));
q=q(reqd);
%-------------------------------------------------------------------------
fprintf('\nNumber of stator turns Ts = ');
disp(round(Ts));
fprintf('\nTotal number of stator slot per phase per pole Ss');
disp(Ss);
fprintf('\nSlot pitch Yss = ');
Yss = pi*D/Ss;
disp(Yss);
Zss = 6*round(Ts);
fprintf('\nTotal Coonductors Zss = ');
disp(Zss);
fprintf('\nNumber of Slots = ');
noofslots = Zss/Ss;
disp(noofslots);
if(choice2==1)
    Es=Vin;
else
    Es=Vin/sqrt(3);
end
Is=Q*1e3/(3*Es);
asc=Is/4;%assuming current density to be 4A/mm2
fprintf('\nCross Sectional Area of Stator Conductor= ');
disp(asc);
fprintf('\n--------------------------------------------------------------');
%%
%Rotor Characteristics calculation
fprintf('\n ROTOR CHARACTERISTICS ');
fprintf('\n --------------------- ');
lg=0.2+2*sqrt(D*L);
fprintf('\nLength of Air Gap lg= ');
disp(lg);
%------------------CALCULATING Sr-------------------------------
no_allow=[0,+p,-p,+2*p,-2*p,+3*p,-3*p,+5*p,-5*p,1,-1,2,-2,p+1,p-1,-p-1,-p+1,p+2,p-2,-p-2,-p+2];
Sr=[Ss:1:1000];
diff=zeros(1,length(Sr));
for i=1:length(Sr)
    diff(i)=Ss-Sr(i);
end %creating the difference array
 
for i=1:1000
    counter=0;
    for j=1:21
        if(diff(i)==no_allow(j))
            counter=counter+1;
        end
    end
    if(counter==0)
        break;
    end
end %calculating the required value of Ss-Sr
Sr=round(Ss+diff(i));
fprintf('\nRotor Slots Sr=');
disp(Sr);
%---------------------Sr CALCULATED-------------------------------
%Diameter of rotor
Dr=D-2*lg;
Lr=L;
fprintf('\nLength of Rotor Lr=');
disp(Lr);
fprintf('\nOuter Diameter of Rotor Dr=');
disp(Dr);
Angle_bw_adj_poles=360/(Ns*p);
fprintf('\nAngle Between Adjacent Poles= ');
disp(Angle_bw_adj_poles);
Electric_angle_skew=(2*pi)/Ns;
fprintf('\nElectrical Angle Skew= ');
disp(Electric_angle_skew);
Ib=(0.85*6*Is*Ts)/Sr; %bar current
fprintf('\nBar Current of Rotor Ib= ');
disp(Ib);
delb=6;%Assumption of 6A/mm2 
ab=Ib/delb;% bar area
fprintf('\nBar Area of Rotor ab=');
disp(ab);
Ie=(Sr*Ib)/(pi*p); %end ring current
dele=7;%Assumption of 7A/mm2
ae=Ie/dele; %end ring area 
ysr=pi*Dr/Sr;
Doer=Dr;
de=10;%depth of the ring is considered 10mm
Dier=Doer-2*de; % inner dia of ring
fprintf('\nInner Diameter of the Rotor= ');
disp(Dier);
Dme=(Doer+Dier)/2; % mean dia
re=0.021*pi*Dme*1e-3/ae;
ohmlossER=2*(Ie^2)*re ;% copper loss in end ring
ohmlossTR=ohmlossER; %neglecting the copper loss at bar
fprintf('\nTotal Copper Loss in rotor= ');
disp(ohmlossTR);
% (rotor copper loss/rotor output)=s/(1-s)
s=ohmlossTR/(Po*1e3+ohmlossTR); % full load slip
fprintf('\nFull Load Slip s=');
disp(abs(s));
fprintf('\n---------------------THE END--------------------');
