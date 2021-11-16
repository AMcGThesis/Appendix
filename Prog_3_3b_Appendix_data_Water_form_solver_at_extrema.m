close all
clear all
%% Wand and nozzle parameters

for jj = 1:3% Loop for wand and nozzle parameters. Here, a 1m scale, 10m scale, and 100m scale wand
    d_loop   =[.00032;.0032;.032];%nozzle jet diameter
    L_loop   =[1;10;100]%wand length
    Vnoz_loop=[2.2;7;22]; %m/s
    
    
    for k = 12:4:28 %Number of nozzles for evaluation
        
        for j = 1:200 %Loop to vary area ratio for a particular scale
            d = d_loop(jj);%2*sqrt(Ax/pi)
            L = L_loop(jj);
            Anoz = (pi/4)*d^2;
            Atotal = k*Anoz
            Apipe = (1+.05*j)*Atotal;
            Ar(j,k) = Apipe / Atotal;
            
            %Initialisation
            D = (4/pi)*sqrt(Apipe);% = pi/4 d^2
            fd = .03; %Darcy Friction Factor
            p = 1000; %kg/m^3
            g = 9.81;
            dh =L/k;
            Lentry = 1.5%m
            n = k+1; %(number of nozzles + 1 (+1 is for the entry length))
            Patm = 1e5; %Pa, absolute
            Qd = .5 %Nozzle Discharge Coefficient
            K = 1/(Qd^2) %Nozzle Resistance coefficient
            
            P = zeros(n,1); %Pa
            Ph = zeros(n,1); %Pa
            pg = zeros(n,1);%Pa
            mpipe = zeros(n,1);%kg/s
            mnoz = zeros(n-1,1); %kg/s
            Vpipe = zeros(n,1); %m/s
            Vnoz = zeros(n-1,1);%m/s
            Pdyn = zeros(n-1,1) %Pa            
            Vnoz(end)=Vnoz_loop(jj); %m/s %sets velocity of the top jet to suit the temporal boundary conditions
            
            %% Initialise flow properties for top nozzle, based on a desired end nozzle velocity
            mnoz(end) = p*Anoz(end)*Vnoz(end);         %End nozzle mass flow rate, based on flow velocity
            
            %Properties in pipe, in order to deliver the requisire amount of water
            %to the top nozzle
            Vpipe(end) = Anoz*Vnoz(end)/Apipe;         %Pipe flow velocity
            mpipe(end) = p*Vpipe(end)*Apipe;           %Pipe mass flow rate
            P(end) = Patm + K*p*Vnoz(end)^2/2;         %Nozzle pressure  = atmospheric + nozzle factor * dynamic pressure at the fastest flow point (the restrictor)
            
            %Pressure delta from top nozzle to next nozzle down the wand
            Ph(end) = dh*fd*Vpipe(end)^2/(2*g*D)*9804; %Pipe Pressure loss per interval(Pascal)
            pg(end) = 0;                               %static pressure gain per nozzle)
            Pdyn(end) = 0;                             %No flow at distal end
            
            %% Evaluate subsequent nozzles, from the flow properties at the top nozzle
            for i = 1:n-2 %for each nozzle other than the top (and leaving room to evaluate an entry length effect)
                P(end-i) =  P(end-i+1)+Pdyn(end-i+1)+Ph(end-i+1)+pg(end-i+1); %P at nozzle i
                
                %P(end) is the pressure AT the top nozzle, P(2) at the bottom nozzle.
                mnoz(end-i) = (Anoz*Qd*p) * sqrt(2*((P(end-i)-Pdyn(end-i))-Patm)/p); %End nozzle mass flow rate at nozzle i
                Vnoz(end-i) = mnoz(end-i)/(p*Anoz);%Velocity of flow through nozzle i

                Vpipe(end-i) = Vnoz(end-i)*Anoz/Apipe + Vpipe(end-i+1); %Pipe flow velocity for section leading up to nozzle i
                mpipe(end-i) = p*Vpipe(end-i)*Apipe; %Pipe mass flow rate fir section leading to nozzle i
 
                Ph(end-i) = dh*fd*Vpipe(end-i+1)^2/(2*g*D)*9804;
                Pdyn(end-i) = (1000/2)*(Vpipe(end-i)^2-Vpipe(end-i+1)^2); %Pumping losses in section leading to nozzle i (pipe + gravity)
                pg(end-i) = p*g*dh;
            end
            PP(j,k) = P(2); %stored vector of max base pressures for evaluating pumping requirements
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
            Re(j,k) = 1000*Vpipe(1)*D / .00131; 
        end
        
        %plot influence of area ratio on pumping pressure requirement for k
        %nozzles
        figure(jj)
        subplot(2,1,1)
        hold on
        plot(Ar(:,k),PP(:,k)/1000, 'linewidth', 1)
        pbaspect([2 1 1])
        xlabel('Area Ratio n*A_p_i_p_e/A_n_o_z','fontsize', 20)
        ylabel('Pump Outlet Pressure (kPa)','fontsize', 20)
        legend('12 Nozzles on Pole','16 Nozzles on Pole', '20 Nozzles on Pole','24 Nozzles on Pole','28 Nozzles on Pole','FontSize',20)
        set(gca, 'fontsize',20)
    end
    
    %plot jet pattern for the case of 20 nozzles.
    figure(jj)
    subplot(2,1,2)
    plot(x',Y','linewidth',1)
    axis([-1.5*L 1.5*L 0 L]);
    pbaspect([2 1 1])
    grid on
    xlabel('Radial Jet Throw Distance (m)')
    ylabel('Height up Wand (m)')
    set(gca, 'fontsize',20)
end