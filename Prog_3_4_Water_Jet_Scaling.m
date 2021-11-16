clc
clear all
close all

%test wand structural properties
l = 0:0.01:15;
g = 9.81; %m/s^2
p = 2020 %kg/m^3
E = 30e9 %GPa
pw = 1000 %kg/m^3
y=3 %gravity parameter

%% Top jet exit velocity
v = l./(sqrt(2*l/g)); %Temporal boundary condition to ensure throw distance = wand length
figure(1)
plot(l,v)
xlabel('Pole Length (m)','FontSize',18)
ylabel(' Required exit velocity of top jet','FontSize',18)
clear l

%% Nozzle jet diameter
do = 3.2;% original size - 3.2mm from water whirler
x1=10; %original wand length from water whirler
x = 0:.01:15; %range of wand lengths (scale ratio for RH = x/1.28)
d = (x / x1)*do; %Nozzle diam over range of scales.
Anoz = (pi/4)*d.^2;

figure(3001)
plot(x/1.28,d)
grid on
xlabel('Scale Ratio (l/l_o)','FontSize',16)
ylabel('Nozzle Diameter[mm]','FontSize',16)


%% Water Critical Diameter vs Minimum bore diameter w a 3:1 area ratio 
i=0;

syms do 
    for l = 1:1:15 %Loop iterates for a range of potential wand lengths
          R =  solve(0 == y*E*do^2+y*E*(do-2*i)^2-16*p*g*l^3); %Crit min diam for an empty wand.                   R = double(R);
          
          %Loop to check roots are suitable, and set minimum bore diameter
          %as crit min. diameter 
          d_crit(1,l) = 0;  
          for k = 1:length(R)
              if abs(imag(R(k))) <= .0001 && R(k)>=0
                   d_crit(1,l) = R(k);
              end      
          end
          L(l) = l;
    end
    
figure(3000)
for i = 1:5 %Loop alters the number of nozzles, 'n'
    n = 8+4*i; 
    Apipe(:,i) = 3*n*Anoz;
    dpipe(:,i) = sqrt((4/pi).*Apipe(:,i));
    hold on
    plot(x/1.28,dpipe(:,i),'-')%'color',c(i,:))
    str(i,:) = sprintf('%d Nozzles on Pole', n);
end    
plot(L/1.28,1000*d_crit, 'k --','linewidth',2)
grid on
xlabel('Scale Ratio (l/l_o)','FontSize',16)
ylabel('Minimum Bore Size [mm]','FontSize',16)

legend(str(1,:),str(2,:),str(3,:),str(4,:),str(5,:),'Critical Diameter','FontSize',16)

