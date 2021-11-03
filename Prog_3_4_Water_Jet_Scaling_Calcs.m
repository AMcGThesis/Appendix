clc
clear all
close all
%%
l = 0:0.01:15;
g = 9.81;
p = 2020
E = 30e9
pw = 1000
y=3
v = l./(sqrt(2*l/g));
figure(1)
plot(l,v)
xlabel('Pole Length (m)','FontSize',18)
ylabel(' Required exit velocity of top jet','FontSize',18)
clear l
%%
%Viewing distance
do = 3.2;%
x1=10;
Aco = (pi/4)*(2*x1*tand(13)).^2;
dco = sqrt(4*Aco/pi);
%Normalise area at 15m to '1' signifying a scaling factor for closer and
%further distances.
x = 0:.01:15;

Ax =(pi/4)*(2*x*tand(13)).^2;
dx = sqrt((4/pi).*Ax);

%d = (dx/dco) * do
d = (x / x1)*do;
%d = sqrt(4*Ax/pi)



Anoz = (pi/4)*d.^2;
%c = jet(5)
figure(3000)
for i = 1:5
    n = 8+4*i;
    Apipe(:,i) = 3*n*Anoz;
    dpipe(:,i) = sqrt((4/pi).*Apipe(:,i));
    hold on
    plot(x/1.28,dpipe(:,i),'-')%'color',c(i,:))
    str(i,:) = sprintf('%d Nozzles on Pole', n);
end
%% Water Critical Diameter
i=0;

syms do 
    for l = 1:1:15
          R =  solve( (16*g*l^3*p)/(E*(do^2 + (do - 2*i)^2)) - y + (16*g*l^3*pw*(do - 2*i)^2)/(E*(do^4 + (do - 2*i)^4)) == 0);
          R2 = solve( (8*g*l^3*p)/(E*do^2) - y + (8*g*l^3*pw)/(E*do^2) == 0);
          R = double(R);
          d_crit(1,l) = 0;  
          for k = 1:length(R)
              if abs(imag(R(k))) <= .0001 && R(k)>=0
                   d_crit(1,l) = R(k);
              end      
          end
          L(l) = l
    end
plot(L/1.28,1000*d_crit, 'k --','linewidth',2)
grid on
xlabel('Scale Ratio (l/l_o)','FontSize',16)
ylabel('Minimum Bore Size [mm]','FontSize',16)

legend(str(1,:),str(2,:),str(3,:),str(4,:),str(5,:),'Critical Diameter','FontSize',16)
figure(3001)
plot(x/1.28,d)
grid on
xlabel('Scale Ratio (l/l_o)','FontSize',16)
ylabel('Nozzle Diameter[mm]','FontSize',16)

