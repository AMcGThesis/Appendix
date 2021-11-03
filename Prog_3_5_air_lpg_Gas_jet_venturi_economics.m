clc
clear all
close all
%d = linspace(.003,.006,.0002)
Ax = (pi/4)*.00435^2;% - (20*(pi/4)*(.0009^2-.0006^2))
%d = 2*sqrt(Ax/pi)
d = linspace(.000,.006,100);
dentry = 2*sqrt( ((pi/4)*.0254^2) / pi );

Aentrytoexit = ((pi/4)*.0254^2)/(20*Ax);
LPG_bottle_cost = 115 ;%Dollars, est. BOC - aug 2018
LPG_bottle_mass = 45; %kg 
LPG_cost_per_kilo= LPG_bottle_cost / LPG_bottle_mass;

for i  = 1:length(d)
[mdot_air(i),mdot_lpg(i),P_air(i),P_LPG(i)] = air_lpg_Isen_nozzle_flow_aug(d(i)) 
end
(25.4^2/.0435^2)
mdot_lpg_one_noz = (mdot_lpg / 20);
mdot_air_one_noz = (mdot_air / 20);
clear mdot_lpg mdot_air
LPG_cost_hour_one_noz = (mdot_lpg_one_noz) * (3600) * LPG_cost_per_kilo;


for i = 12:4:28
    mdot_lpg(i,:) = 3600*mdot_lpg_one_noz*i;
    mdot_air(i,:) = 3600*mdot_air_one_noz*i;
    LPG_cost_hour(i,:) = LPG_cost_hour_one_noz*i;
    
figure(1)
 hold on
    plot(d,mdot_lpg(i,:),'linewidth',1);%LPG_cost_hour(i,:))
    xlabel('Venturi Diameter [mm]', 'fontsize', 18)
    ylabel('LPG Flow Rate (kg/hr at 50kPa)', 'fontsize', 18)
end

axis([0 6e-3 0 mdot_lpg(i,end)])
    yyaxis right
    %plot(d,mdot_lpg*3600)
    ylabel('LPG Cost / Hour ($)', 'fontsize', 18)
   ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
    %legend('LPG (50kPa)','Flow rate kg/hr')
    axis([0 6e-3 0 LPG_cost_hour(i,end)])
    legend('12 Nozzles', '16 Nozzles', '20 Nozzles', '24 Nozzles', '28 Nozzles', 'fontsize', 14)
    grid on