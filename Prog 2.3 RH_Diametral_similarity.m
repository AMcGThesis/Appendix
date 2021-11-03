%Similar diametral pairs for Rotating Harmonic
clc
clear all
close all
%% Geometric and Structural property definitions
%Properties of scaled wand
l = 0:.1:15; %Test lengths for scaling
y=1.6; %Gravity parameter of Rotating Harmonic
E=30e9; %Pa Young's mod. (GFRP test wand)
p = 2020; %kg/m^3 Desity (GFRP)

g = 9.81; %ms^2

%properties of original wand
doo = .0028; % m OD of RH's steel wand
doi = 0; % m ID of RH's steel wand
po = 7850; %kg/m^3 (RH is steel)
lo = 1.28; %m Wand length (RH)
Sig2=691e6; %Pa Stress from RH STRUCTURAL LOADS CALC
So = 5.44   %N Shear from RH  STRUCTURAL LOADS CALC
Ao = .175 %m RH operating amplitude. 
f2 = 7.64 %2nd mode resfreq of rotating harmonic

%% 
SigP = [1,1.8/1.38,2.5/1.38,3/1.38]; %Test stress conc. values
tvec = [5,10,20]*10^-3; % mm test wall thicknesses


%% Diameters
figure(2)
syms di
hold on
d = sqrt(16*p*g*l.^3 / (2*E*y)); %m crit. minimum diameter (eqn 2.3.20)
for j = 0:1:length(tvec)-1
    i = tvec(j+1);
    for l = 1:1:15
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


l = [1:1:15]; %Test lengths for scaling
    hold on
    plot(l(2:end)/lo, 1000*Do(1,2:end))
    plot(l(3:end)/lo,1000*Do(2,3:end))
    plot(l(5:end)/lo,1000*Do(3,5:end))

    set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Outer Diameter [mm]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')

%% Stresses
figure(3)
for j = 1:length(tvec)
    for k = 1:length(l)
        Sig(j,k) = (p*l(k)^2*Do(j,k)*(doo^2+doi^2))/(po*lo^2*doo*(Do(j,k)^2+Di(j,k)^2)); %Stress ratio: (2.3.7)
    end
end
Sig = Sig/10^6 %Conversion into MPa
hold on
   plot(l(:)/lo, Sig2*Sig(1,:))
    plot(l(2:end)/lo,Sig2*Sig(2,2:end))
    plot(l(4:end)/lo,Sig2*Sig(3,4:end))

    set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Stress [MPa]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')
%% Determine max amplitude supported by 1/3rd yield criterion
figure (303)

l = [1:1:15];
hold on
for k = 1:length(SigP) %Eval. supported amplitudefor each stress conc (k=1 is base stress)
    for j = 1:length(tvec) %for each test wall thickness
        for i = 1:length(l) %for each scale ratio
            Amp{k}(j,i) = (Ao/lo)*267/(SigP(k)*Sig2*Sig(j,i)) 
            %Dimnless amplitude required for 1/3rd yield stress. 
            %(800MPa is ~ yield for a pultrusion w wrapped material. 
        end
    end
end

    plot(l(3:end)/lo,   Amp{1}(1,3:end)    ,'linewidth', 1.5,'color', [0, 0.4470, 0.7410] )
    plot(l(5:end)/lo,   Amp{1}(2,5:end)    ,'linewidth', 1.5,'color', [0.8500, 0.3250, 0.0980] )
    plot(l(7:end)/lo,   Amp{1}(3,7:end)    ,'linewidth', 1.5,'color', [0.9290, 0.6940, 0.1250] )
    plot(linspace(0,12,length(l)),0.14*ones(length(l),1),'k --','linewidth', 1.5)
 
    plot(l(3:end)/lo,   Amp{3}(1,3:end)    ,'linewidth', 1.5,'color', [0, 0.4470, 0.7410] )
    plot(l(5:end)/lo,   Amp{3}(2,5:end)    ,'linewidth', 1.5,'color', [0.8500, 0.3250, 0.0980] )
    plot(l(7:end)/lo,   Amp{3}(3,7:end)    ,'linewidth', 1.5,'color', [0.9290, 0.6940, 0.1250] )%   
    
    set(gca,'FontSize',24)
    xlabel('Scale Ratio [l/l_o]')
    ylabel('Dimensionless Amplitude [m/m]')
    legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness','Geometrically Similar Amplitude')

figure(301)

hold on
    plot(l(:)    /lo,   SigP(1)*Sig2*Sig(1,:)    ,'color', [0, 0.4470, 0.7410]   )
    plot(l(2:end)/lo,   SigP(1)*Sig2*Sig(2,2:end),'color', [0.8500, 0.3250, 0.0980] )
    plot(l(4:end)/lo,   SigP(1)*Sig2*Sig(3,4:end),'color', [0.9290, 0.6940, 0.1250] )
 
    plot(l(:)    /lo,   SigP(2)*Sig2*Sig(1,:)   ,'--','color', [0, 0.4470, 0.7410]   )
    plot(l(2:end)/lo,   SigP(2)*Sig2*Sig(2,2:end),'--','color', [0.8500, 0.3250, 0.0980] )
    plot(l(4:end)/lo,   SigP(2)*Sig2*Sig(3,4:end),'--','color', [0.9290, 0.6940, 0.1250] )
    
    plot(l(:)    /lo,   SigP(3)*Sig2*Sig(1,:)    ,'-','color', [0, 0.4470, 0.7410]   )
    plot(l(2:end)/lo,   SigP(3)*Sig2*Sig(2,2:end),'-','color', [0.8500, 0.3250, 0.0980] )
    plot(l(4:end)/lo,   SigP(3)*Sig2*Sig(3,4:end),'-','color', [0.9290, 0.6940, 0.1250] )
    
    plot(l(:)    /lo,   SigP(4)*Sig2*Sig(1,:)    ,'--','color', [0, 0.4470, 0.7410]   )
    plot(l(2:end)/lo,   SigP(4)*Sig2*Sig(2,2:end),'--','color', [0.8500, 0.3250, 0.0980] )
    plot(l(4:end)/lo,   SigP(4)*Sig2*Sig(3,4:end),'--','color', [0.9290, 0.6940, 0.1250] )
     set(gca,'FontSize',20)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Stress [MPa]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')


%% Moments
Mo = Sig2*(pi/64)*(doo^4)/.0014% Bending Stress = My/I 

figure(4)
for j = 1:length(tvec)
    for k = 1:length(l)
       M(j,k) = (p*(Do(j,k)^2-Di(j,k)^2)*l(k)^2)/(po*(doo^2-doi^2)*lo^2); %Moment ratio (2.3.5)
    end
end
M = M/10^3; %Convert to kNm
hold on
    plot(l(:)/lo, Mo*M(1,:))
    plot(l(2:end)/lo,Mo*M(2,2:end))
    plot(l(3:end)/lo,Mo*M(3,3:end))
 
    set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Moment [kNm]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')

%% Shear

l = [1:1:15];
figure(5)
for j = 1:length(tvec)
    for k = 1:length(l)
        S(j,k) = (p*(Do(j,k)^2-Di(j,k)^2)*l(k))/(po*(doo^2-doi^2)*lo); %Shear ratio (2.3.10)
    end
end
S = S/10^3; %Convert to kN
hold on
    plot(l(:)/lo, So*S(1,:))
    plot(l(3:end)/lo,So*S(2,3:end))
    plot(l(5:end)/lo,So*S(3,5:end))
    set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Shear Force [kN]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')


%% Frequency

figure(6)
l = linspace(1,15,1000);
for j = 1:length(l)
    f(j) = sqrt(lo/l(j)); %Frequency ratio (2.3.13)
end
    plot(l(:)/lo, f2*f(1,:))
    
    set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Frequency [Hz]')
legend('2nd Bending Mode')


