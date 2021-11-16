%Similar diametral pairs for Zebra
clc
clear all
close all
%% Geometric and Structural property definitions
%Properties of scaled wand
l = 0:.1:8*2.1; %Test wand lengths
y=4.59;      %Gravity parameter of Zebra
E=30e9;      %Pa Young's mod. (GFRP test wand)
p = 2020;    %kg/m^3 Desity (GFRP)
g = 9.81;    %ms^2

%Properties of original wand
doo = .004;  % m OD of Zebras OG GFRP wand
doi = 0;     % m ID of RH's steel wand
po = 2020;   %kg/m^3 (RH is steel)
lo = 2.1;    %m Wand length (RH)
Sig2 = 139.3;%Pa Stress from Zebra (2nd mode)STRUCTURAL LOADS CALC
Sig3 = 362.9;% " (3rd mode)
M2 = .875;   %Nm Moment from Zebra (2nd mode)STRUCTURAL LOADS CALC
M3 = 2.28;   % " (3rd mode)
S2 = 1.87    %N Shear from Zebra   (2nd mode)STRUCTURAL LOADS CALC
S3 = 8.29    % " (3rd mode)
f2 = 3.39    %Hz res. freq Zebra   (2nd mode) STRUCTURAL LOADS CALC
f3 = 9.76    % " (3rd mode)
%% initialise
 tvec = 2*[2.5,5,10]*10^-3; %Test wall thicknesses (double)
 SigP = [1,1.8,2.5,3]; %Wand SCF's
%% Diameters
figure(2)
syms di
hold on
d = sqrt(16*p*g*l.^3 / (2*E*y)); %m crit. minimum diameter (eqn 2.3.20)

for j = 0:1:length(tvec)-1
    i = tvec(j+1);
    for l = 1:1:16
      R =  solve( y*E*(di+i)^2+y*E*di^2 == 16*p*g*l^3); %Roots of eqn 2.3.19)
      R = double(R);
      if imag(R(2)) <= .0001 && R(2)>=0 %Find diametrally similar pairs (real+positive only)
            Di(j+1,l) = R(2);
            Do(j+1,l) = Di(j+1,l)+i;
      else
            Di(j+1,l) = 0;
            Do(j+1,1) = 0;      
      end
    end
     t(j+1) = i;
end

l = [1:1:16];
    hold on
    plot(l(3:end)/lo, 1000*Do(1,3:end))
    plot(l(5:end)/lo,1000*Do(2,5:end))
    plot(l(7:end)/lo,1000*Do(3,7:end))

set(gca,'FontSize',14)
grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Outer Diameter [mm]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')

%% Stresses
figure(3)
for j = 1:length(tvec)
    for k = 1:length(l)
           Sig(j,k) = (p*l(k)^2*Do(j,k)*(doo^2+doi^2))/(po*lo^2*doo*(Do(j,k)^2+Di(j,k)^2));%Stress ratio: (2.3.7)
    end
end

hold on
    plot(l(3:end)/lo, Sig3*Sig(1,3:end),'-','color', [0, 0.4470, 0.7410] )
    plot(l(5:end)/lo, Sig3*Sig(2,5:end),'-','color', [0.8500, 0.3250, 0.0980]	 )
    plot(l(7:end)/lo, Sig3*Sig(3,7:end),'-','color', [0.9290, 0.6940, 0.1250])
    
set(gca,'FontSize',14)
grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Third Mode Stress [MPa]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')

%% Determine max amplitude supported by 1/3rd yield criterion
l = [1:1:16];
figure (303)
hold on
for k = 1:length(SigP) %Eval. supported amplitude for each stress conc (k=1 is base stress)
    for j = 1:length(tvec) %for each test wall thickness
        for i = 1:length(l) %for each scale ratio
            Amp{k}(j,i) = (.35/lo)*267/(SigP(k)*Sig3*Sig(j,i)) 
         %Dimnless amplitude required for 1/3rd yield stress. 
        end
    end
end

    plot(l(3:end)/lo,   Amp{1}(1,3:end)    ,'linewidth', 1.5,'color', [0, 0.4470, 0.7410] )
    plot(l(5:end)/lo,   Amp{1}(2,5:end)    ,'linewidth', 1.5,'color', [0.8500, 0.3250, 0.0980] )
    plot(l(7:end)/lo,   Amp{1}(3,7:end)    ,'linewidth', 1.5,'color', [0.9290, 0.6940, 0.1250] )

    plot(l(3:end)/lo,   Amp{2}(1,3:end)    ,'linewidth', 1.5,'color', [0, 0.4470, 0.7410] )
    plot(l(5:end)/lo,   Amp{2}(2,5:end)    ,'linewidth', 1.5,'color', [0.8500, 0.3250, 0.0980] )
    plot(l(7:end)/lo,   Amp{2}(3,7:end)    ,'linewidth', 1.5,'color', [0.9290, 0.6940, 0.1250] )%    
    
    plot(l(3:end)/lo,   Amp{3}(1,3:end)    ,'linewidth', 1.5,'color', [0, 0.4470, 0.7410] )
    plot(l(5:end)/lo,   Amp{3}(2,5:end)    ,'linewidth', 1.5,'color', [0.8500, 0.3250, 0.0980] )
    plot(l(7:end)/lo,   Amp{3}(3,7:end)    ,'linewidth', 1.5,'color', [0.9290, 0.6940, 0.1250] )%    
    
    plot([0 100], [.17,.17],'k --','linewidth', 1.5)
set(gca,'FontSize',28)
xlabel('Scale Ratio [l/l_o]')
ylabel('Dimensionless Amplitude [m/m]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness', 'Geometrically Similar Amplitude')
ylim([0.02 0.18])
xlim([0 8])
grid minor
%% Moments
figure(4)
for j = 1:length(tvec)
    for k = 1:length(l)
        M(j,k) = (p*(Do(j,k)^2-Di(j,k)^2)*l(k)^2)/(po*(doo^2-doi^2)*lo^2); %Moment ratio (2.3.5)
    end
end
M = M/10^3;
hold on
    plot(l(:)    /lo,M3*M(1,:),'-','color', [0, 0.4470, 0.7410] )
    plot(l(2:end)/lo,M3*M(2,2:end),'-','color', [0.8500, 0.3250, 0.0980]	 )
    plot(l(3:end)/lo,M3*M(3,3:end),'-','color', [0.9290, 0.6940, 0.1250])
     set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Moment [kNm]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')

%% Shear
l = [1:1:16];
figure(5)
for j = 1:length(tvec)
    for k = 1:length(l)
         S(j,k) = (p*(Do(j,k)^2-Di(j,k)^2)*l(k))/(po*(doo^2-doi^2)*lo); %Shear ratio (2.3.10)
    end
end
S = S/10^3;
hold on
 
    plot(l(:)    /lo,S3*S(1,:),'-','color', [0, 0.4470, 0.7410] )
    plot(l(4:end)/lo,S3*S(2,4:end),'-','color', [0.8500, 0.3250, 0.0980]	 )
    plot(l(7:end)/lo,S3*S(3,7:end),'-','color', [0.9290, 0.6940, 0.1250])
    set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Shear Force [kN]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')


%% Frequency
figure(6)
l = linspace(1,16,1000);
for j = 1:length(l)
    
f(j) = sqrt(lo/l(j));  %Frequency ratio (2.3.13)
end
    plot(l(:)/lo, f2*f(1,:),l(:)/lo, f3*f(1,:)) 
    set(gca,'FontSize',14)
    grid on
    axis([0 8 0 15.1])
xlabel('Scale Ratio [l/l_o]')
ylabel('Frequency [Hz]')
legend('2nd Bending Mode', '3rd Bending Mode')

