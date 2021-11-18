%Evaluates scaled properties using the outputs of Prog 3_2.
clc
clear all
close all
%% Diameters
figure(2)
clear l
syms E y di do pw i l p g 

%Using 'test wand' properties
l = 0:.1:20;
y=3; %Desired gravity parameter 
E=30e9; %GPa Young's mod. (GFRP test wand)
p = 2020; %kg/m^3 Desity (GFRP)
g = 9.81; %ms^2
pw = 1000 %kg/m^3 water density
syms di

%Unit wand properties from prog. 3.1 to scale with respect to
doo = 7.23%.0028;
doi = 0;
po = 1%7850;
lo = 1%1.28;
Mo = 517%
So = 2369; %Second mode shear of FH unit canti
f2 = 6.17 %2nd mode moment of FH unit canti

tvec = [2.5,5,10]*10^-3; %Wall thickness range to be investigated

%% Statically similar diameters
for j = 0:1:length(tvec)-1
    i = tvec(j+1);
      
    for l = 1:1:15     
      R =  solve( (16*g*l^3*p)/(E*(do^2 + (do - 2*i)^2)) - y + (16*g*l^3*pw*(do - 2*i)^2)/(E*(do^4 - (do - 2*i)^4)) == 0); %Static similarity criterion for fluid filled wands
      R = double(R)
      Di(j+1,l) = 0;
      Do(j+1,1) = 0;  
      for k = 1:length(R)
          if abs(imag(R(k))) <= .0001 && R(k)>=0%Select only real, positive rools
               Do(j+1,l) = R(k);
               Di(j+1,l) = Do(j+1,l)-2*i;
          end      
      end
    end
     t(j+1) = i;
end

% plotting potential wand size choices
l = [1:1:15];
    hold on
    plot(l(3:end)/1.28, 1000*Do(1,3:end))
    plot(l(4:end)/1.28,1000*Do(2,4:end))
    plot(l(6:end)/1.28,1000*Do(3,6:end))
    set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Water Filled Outer Diameter [mm]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')
axis([0 12 0 200])

l = [1:1:15];

%% Grav param. change when draining the Y=3 wand
figure(10)
for j = 1:length(tvec)
    for k = 1:length(l)
       y_test(j,k) = 16*p*(Do(j,k)^2-Di(j,k)^2)*g*l(k)^3 / (E*(Do(j,k)^4-Di(j,k)^4));
    end
end

hold on
    plot(l(3:end)/1.28, y_test(1,3:end),'linewidth',1.25)
    plot(l(4:end)/1.28,y_test(2,4:end),'linewidth',1.25)
    plot(l(6:end)/1.28,y_test(3,6:end),'linewidth',1.25)
    set(gca,'FontSize',18)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Gravity Parameter of a Gas Performance y_g')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness')


%% Moments Water
figure(8)
for j = 1:length(tvec)
    for k = 1:length(l)
       M(j,k) =  ((pw*Di(j,k)^2+ p*(Do(j,k)^2-Di(j,k)^2))*l(k)^2) / ...
                 ((pw*doi^2    + po*(doo^2-doi^2))        *lo^2); %Moment ratio
    end
end
M = Mo*M;
hold on
    plot(l(3:end)/1.28,M(1,3:end)/10^3)
    plot(l(4:end)/1.28,M(2,4:end)/10^3)
    plot(l(6:end)/1.28,M(3,6:end)/10^3)
    set(gca,'FontSize',14)
    grid on
     axis([0 12 0 75])
xlabel('Scale Ratio [l/l_o]')
ylabel('Water Moment [kNm]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness','25mm wall thickness')

%% Stress Water
figure(7)
for j = 1:length(tvec)
    for k = 1:length(l)
        SigW(j,k) = M(j,k)*Do(j,k) / (2*(pi/64)*(Do(j,k)^4-Di(j,k)^4));
      end
end
hold on
    plot(l(3:end)/1.28,SigW(1,3:end)/10^6)
    plot(l(4:end)/1.28,SigW(2,4:end)/10^6)
    plot(l(6:end)/1.28,SigW(3,6:end)/10^6)
    set(gca,'FontSize',14)
    grid on
    axis([0 12 0 800])
xlabel('Scale Ratio [l/l_o]')
ylabel('Water Filled Stress [MPa]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness','25mm wall thickness')

%% Shear Water
figure(9)
for j = 1:length(tvec)
    for k = 1:length(l)
      S(j,k) =  (pw*Di(j,k)^2+p*(Do(j,k)^2-Di(j,k)^2))*l(k) / ...
              ((pw*doi^2    +po*(doo^2-doi^2))         *lo);
    end
end
S = S/10^3;
hold on
    plot(l(3:end)/1.28,So*S(1,3:end))
    plot(l(4:end)/1.28,So*S(2,4:end))
    plot(l(6:end)/1.28,So*S(3,6:end))
    set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Water Shear Force [kN]')
legend('2.5mm wall thickness','5mm wall thickness','10mm wall thickness','25mm wall thickness')
axis([0 12 0 20])


%% Frequency Water
figure(6)
l = linspace(1,15,1000);
for j = 1:length(l)
    
f(j) = sqrt((lo/l(j)));
end
    plot(l(:)/1.28, f2*f(1,:))
    
    set(gca,'FontSize',14)
    grid on
xlabel('Scale Ratio [l/l_o]')
ylabel('Water Frequency [Hz]')
legend('2nd Bending Mode')

