%Water specifications used to select a pump in Chapter 6.
close all
clear all
%% Wand and nozzle parameters
d   = .0032; %nozzle jet diameter
L   =[10]    %wand length
Vtopjet=[7]; %m/s
k = 20       %Th number of nozzles to be used
Qd = .5 %Nozzle Discharge Coefficient

%Initialisation. 
Anoz = (pi/4)*d^2;
Atotal = k*Anoz
Apipe = (pi/4)*(.0276^2)%
D = (4/pi)*sqrt(Apipe);% (A= (pi/4) d^2)
fd = .03; %Darcy Friction Factor
p = 1000; % water densitykg/m^3
g = 9.81;
dh =L/k;
Lentry = 1.5%m
n = k+1; %(number of nozzles + 1 (+1 is for the entry length))
Patm = 1e5; %Pa, absolute atmospheric pressure
K = 1/(Qd^2) %Nozzle Resistance coefficient

% 
P = zeros(n,1); %Pa
Ph = zeros(n,1); %Pa
pg = zeros(n,1);%Pa
mpipe = zeros(n,1);%kg/s
mnoz = zeros(n-1,1); %kg/s
Vpipe = zeros(n,1); %m/s
Vnoz = zeros(n-1,1);%m/s
Pdyn = zeros(n-1,1); %Pa
Vnoz(end)=Vtopjet; %m/s %sets velocity of the top jet to suit the temporal boundary conditions

%% Initialise flow properties for top nozzle, based on a desired end nozzle velocity
mnoz(end) = p*Anoz(end)*Vnoz(end); %End nozzle mass flow rate, based on flow velocity

%Properties in pipe, in order to deliver the requisire amount of water
%to the top nozzle
Vpipe(end) = Anoz*Vnoz(end)/Apipe; %Pipe flow velocity
mpipe(end) = p*Vpipe(end)*Apipe;   %Pipe mass flow rate
P(end) = Patm + K*p*Vnoz(end)^2/2; %Nozzle pressure  = atmospheric + nozzle factor * dynamic pressure at the fastest flow point (the restrictor)

%Pressure delta from top nozzle to next nozzle down the wand
Ph(end) = dh*fd*Vpipe(end)^2/(2*g*D)*9804; %Pipe Pressure loss per interval(Pascal)
pg(end) = 0;                               %static pressure gain per nozzle)
Pdyn(end) = 0;                             %No flow at distal end

%% Evaluate subsequent nozzles, from the flow properties at the top nozzle
for i = 1:n-2 %for each nozzle other than the top (and leaving room to evaluate an entry length effect)
    P(end-i) =  P(end-i+1)+Pdyn(end-i+1)+Ph(end-i+1)+pg(end-i+1); % P at nozzle i-1
    
    mnoz(end-i) = (Anoz*Qd*p) * sqrt(2*((P(end-i)-Pdyn(end-i))-Patm)/p); %End nozzle mass flow rate at nozzle i-1
    Vnoz(end-i) = mnoz(end-i)/(p*Anoz);%Velocity of flow through nozzle i-1
    
    Vpipe(end-i) = Vnoz(end-i)*Anoz/Apipe + Vpipe(end-i+1); %Pipe flow velocity for section leading up to nozzle i-1
    mpipe(end-i) = p*Vpipe(end-i)*Apipe; %Pipe mass flow rate fir section leading to nozzle i-1
    
    %Pumping losses in section leading to nozzle i (pipe + gravity + dynamic pressure gain)
    Ph(end-i) = dh*fd*Vpipe(end-i+1)^2/(2*g*D)*9804;
    Pdyn(end-i) = (1000/2)*(Vpipe(end-i)^2-Vpipe(end-i+1)^2);
    pg(end-i) = p*g*dh;
    
end

Vpipe(1) = Vnoz(end-i)*Anoz/Apipe + Vpipe(end-i+1); %Pipe flow velocity for section leading up to nozzle i


%% Nozzle kinematics
t=zeros(1,k);
y=zeros(1,k);
wand_height = L; %metres
dh;
if k == 20
    for i = 1:k %njets=  number of jets
        y(i)=i*dh;
        t(i) = sqrt(2*y(i)/g);  %flight time for a jet i*jetspacing metres above ground
        T(i,:)=linspace(0,t(i),100);
        Y(i,:)=-g*(T(i,:).^2)/2+ y(i); %Y coordinates of water jets)
        x(i,:) = Vnoz(i)*T(i,:); %x-coords of water jets
    end
    
    for i = 1:n
        if mod(i,2) == 0
            x(i,:)=-x(i,:);
            
        end
    end
end
t_bottom = t(1);
t_top =  t(end);


Ar(k) = Apipe / Atotal; %area ratio
Re(k) = 1000*Vpipe(1)*D / .00131; %Reynolds at base. NB .00131 = dynamic viscousity for water ~20 degrees.

%Figure of flow form.
figure(1)
plot(x',Y','linewidth',1)
axis([-1.5*L 1.5*L 0 L]);
pbaspect([2 1 1])
grid on
xlabel('Radial Jet Throw Distance (m)')
ylabel('Height up Wand (m)')
set(gca, 'fontsize',20)

%% Output loop
fprintf('Flow Outputs\n***************************');

fprintf('\n Pressure at base is %4.0f kPa',P(2)/1000)
fprintf('\n Mass flow rate at base is %4.1f kg/s',mpipe(2))

fprintf('\n ***************************');
