%% Satellite parameters
sat.Is = [0.103682,-0.001094, 0.015519;...
         -0.001094, 0.104158,-2.57e-6;...
          0.015519,  2.57e-6, 0.088588];
sat.mass = 40; %Kg

%% Actuators parameters
actuators.magnetorquer.power            = 1536E-3;       %W
actuators.magnetorquer.MaxPower         = 13200-3;       %W
actuators.magnetorquer.voltage          = 12;            %V
actuators.magnetorquer.dimensions       = [250,67,35];   %mm
actuators.magnetorquer.nominalDipole    = 15;            %Am^2
actuators.magnetorquer.maxNominalDipole = 30;            %Am^2
actuators.magnetorquer.A  = actuators.magnetorquer.dimensions(2)*actuators.magnetorquer.dimensions(3)*1E-6;   %Transversal Area (m ^2)
actuators.magnetorquer.n  = actuators.magnetorquer.nominalDipole*actuators.magnetorquer.voltage/(actuators.magnetorquer.A*actuators.magnetorquer.power);           %Number of turns
actuators.magnetorquer.maxCurrent = actuators.magnetorquer.maxNominalDipole/actuators.magnetorquer.voltage;  %Max current (mA)

%% Sensor parameters
sensors.mag.desvEst = 50.00e-9; %T
sensors.mag.res     = 31.25e-9; %T
sensors.mag.normDesvEst =  0.5;

%% Planet parameters
%%%%This will define all of our planet parameters
earth.Radius       = 6378E3;               %% Planet Radius (meters)
earth.Mass         = 5.972E24;             %% Planet Mass (kg)
earth.GravityConst = 6.674E-11;            %% Planet gravity constan (Nm^2/kg^2)
earth.mu           = earth.GravityConst*earth.Mass;        %% G*M product (Nm^2/kg)

%% Orbit Parameters
%%%Orbital parameters
orbit.altitude    = 500E3;   
orbit.semiMajorAxis = orbit.altitude+earth.Radius;
orbit.eccentricity = 0;
orbit.inclination = 97.34; 
orbit.rightAscensionOfAscendingNode = 313.131; 
orbit.argumentOfPeriapsis = 0; 
orbit.trueAnomaly = 0;
orbit.period      = 2*pi/sqrt(earth.mu)*orbit.semiMajorAxis^(3/2); % Orbital period
orbit.vcircular   = sqrt(earth.mu/orbit.semiMajorAxis);            % Circular velocity

%% External disturbances
disturbance = @(t) simplifiedDisturbances(t);
kfunConstant = @(dip)constantK(dip);
kfunVariable = @(dip)variableK(dip);

%% Initial Conditions
%%% Intitial Conditions for Attitude and Angular Velocity (Euler angles)
initial.atitude.rpy0_deg(1)    = 0;           %  Initial Roll  (deg)
initial.atitude.rpy0_deg(2)    = 0;           %  Initial Pitch (deg)
initial.atitude.rpy0_deg(3)    = 0;           %  Initial Yaw   (deg)
%%% Initial Euler Angles and quaternions
initial.atitude.rpy0_rad    = deg2rad([initial.atitude.rpy0_deg(1),...
    initial.atitude.rpy0_deg(2),...
    initial.atitude.rpy0_deg(3)]');       % Euler Angles (rad)
initial.atitude.q0123_0 = EulerAngles2Quaternions(initial.atitude.rpy0_rad);  % Quaternions
%%% Initial angular rates in the body frame (rad/s)
initial.omega.omega0_x = -1*deg2rad(3);
initial.omega.omega0_y = deg2rad(4);
initial.omega.omega0_z = -1*deg2rad(5);

%%Magnetic Field
settings.B_earth_field.model_epoch = '2020';
settings.B_earth_field.decimal_year = 2020 ;
settings.wgs84 = wgs84Ellipsoid;

%% Simulation Parameters
%%% Setting Time Window (just 5 orbit)
settings.startTime = datetime(2023,08,30,12,35,38);
settings.sampleTime = 100;
settings.number_of_orbits = 5;
settings.tfinal = orbit.period*settings.number_of_orbits;
%settings.tfinal = 2;
settings.stopTime = settings.startTime  + seconds(settings.tfinal);

%%% Satellite initial states
settings.X0 = [initial.atitude.q0123_0; initial.omega.omega0_x; ...
               initial.omega.omega0_y ; initial.omega.omega0_z];

%% Create satellite scenario to use SPG4 propagator
sc = satelliteScenario(settings.startTime ,settings.stopTime,settings.sampleTime);
sat.satSGP4 = satellite(sc, orbit.semiMajorAxis, orbit.eccentricity, orbit.inclination, ...
        orbit.rightAscensionOfAscendingNode, orbit.argumentOfPeriapsis, orbit.trueAnomaly);
 
%% Thirty simulation for the k constant
num = 30;

simulationc = struct(); simulationsC = struct();
simulationV = struct(); simulationsV = struct();

%% Create initial conditions
X0 = zeros(7,num);
for k = 1:num
    X0(:,k) =[initial.atitude.q0123_0; deg2rad((rand(1) * 10)-5); ...
                  deg2rad((rand(1) * 10)-5); deg2rad((rand(1) * 10)-5)]; 
end

%% Simulations for k constant
for k = 1:num
    simulationc = detumblingSimulation(simulationc,kfunConstant, ...
                       sat, disturbance, sensors, settings, actuators, X0(:,k));
    simulationc.X0 = X0(:,k);
    if k == 1
        simulationsC = simulationc;
        figure()
    else
        simulationsC = [simulationsC;simulationc];
    end
    %% Plots
    % Rename index
    q0123out = simulationc.x(:,1:4);                           %%Satellite quaternions
    ptpout = Quaternions2EulerAngles(q0123out);    %%Satellite euler angles
    pqrout = simulationc.x(:,5:7);                             %%Angular Rates
    tout = simulationc.tout;
p1 = subplot(2,1,1);
    h1 = plot(tout/3600,rad2deg(pqrout(:,1)),'b-','LineWidth',2); hold on; grid on
    h2 = plot(tout/3600,rad2deg(pqrout(:,2)),'r-','LineWidth',2);
    h3 = plot(tout/3600,rad2deg(pqrout(:,3)),'g-','LineWidth',2);
    max = 0.1*ones(1,length(tout));
    h4 = plot(tout/3600,max,'k--','LineWidth',1);
    h5 = plot(tout/3600,-1*max,'k--','LineWidth',1);
    xlabel('Time (h)'); ylabel('\omega (deg/s)');
    legend([h1,h2,h3,h4,h5],'\omega_x','\omega_y','\omega_z','Tolerance','Tolerance','Location','eastoutside');
    title('Satellite angular rate (deg/s)');

p2 = subplot(2,1,2);
    h1 = plot(tout/3600,rad2deg(simulationc.magnOmega),'LineWidth',2); hold on; grid on;
    max = 3*0.1*ones(1,length(tout));
    h2 = plot(tout/3600,max,'k--','LineWidth',2);
    xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
    title('Satellite angular rate magnitude (deg/s)');
    %xlim([0.0 5.0]);
    legend([h1,h2],'|\omega|','Tolerance','Location','eastoutside');
    linkaxes([p1,p2],'x')
end

%% Simulations for k variable
for k = 1:num
    simulationV = detumblingSimulation(simulationV,kfunVariable, ...
                       sat, disturbance, sensors, settings, actuators, X0(:,k));
    simulationV.X0 = X0(:,k);
    if k == 1
        simulationsV = simulationV;
        figure()
    else
        simulationsV = [simulationsV;simulationV];
    end
    %% Plots
    % Rename index
    q0123out = simulationV.x(:,1:4);                           %%Satellite quaternions
    ptpout = Quaternions2EulerAngles(q0123out);    %%Satellite euler angles
    pqrout = simulationV.x(:,5:7);                             %%Angular Rates
    tout = simulationV.tout;
p1 = subplot(2,1,1);
    h1 = plot(tout/3600,rad2deg(pqrout(:,1)),'b-','LineWidth',2); hold on; grid on
    h2 = plot(tout/3600,rad2deg(pqrout(:,2)),'r-','LineWidth',2);
    h3 = plot(tout/3600,rad2deg(pqrout(:,3)),'g-','LineWidth',2);
    max = 0.1*ones(1,length(tout));
    h4 = plot(tout/3600,max,'k--','LineWidth',1);
    h5 = plot(tout/3600,-1*max,'k--','LineWidth',1);
    xlabel('Time (h)'); ylabel('\omega (deg/s)');
    legend([h1,h2,h3,h4,h5],'\omega_x','\omega_y','\omega_z','Tolerance','Tolerance','Location','eastoutside');
    title('Satellite angular rate (deg/s)');
    %xlim([0.0 5.0]);

p2 = subplot(2,1,2);
    h1 = plot(tout/3600,rad2deg(simulationV.magnOmega),'LineWidth',2); hold on; grid on;
    max = 3*0.1*ones(1,length(tout));
    h2 = plot(tout/3600,max,'k--','LineWidth',2);
    xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
    title('Satellite angular rate magnitude (deg/s)');
    %xlim([0.0 5.0]);
    legend([h1,h2],'|\omega|','Tolerance','Location','eastoutside');
    linkaxes([p1,p2],'x')
end

%% Plots k constant
figure()
%%%Set toletances
max_tolerance      = 0.1*ones(1,length(simulationsC(1).tout));
max_tolerance_norm = 0.3*ones(1,length(simulationsC(1).tout));
%%%Calculate average settlement time 
settlementTime = 0;
for k = 1:num
    %%% Sum SettlementTime
    settlementTime = settlementTime + simulationsC(k).tiempo_inicio;
    %%% Rename index
    pqrout = simulationsC(k).x(:,5:7);      %% Angular Rates
    tout = simulationsC(k).tout;            %% time obtained
    
    if k == 1
        p1 = subplot(2,1,1);
            h1 = plot(tout/3600,rad2deg(pqrout(:,1)),'b-','LineWidth',2); hold on; grid on
            h2 = plot(tout/3600,rad2deg(pqrout(:,2)),'r-','LineWidth',2);
            h3 = plot(tout/3600,rad2deg(pqrout(:,3)),'g--','LineWidth',2);
            h4 = plot(tout/3600,max_tolerance,'k--','LineWidth',1);
            h5 = plot(tout/3600,-max_tolerance,'k--','LineWidth',1);
            xlabel('Time (h)'); ylabel('\omega (deg/s)');
            title('Simulations k constant: Satellite angular rate (deg/s)');
            legend([h1,h2,h3,h4,h5],'\omega_x','\omega_y','\omega_z','Tolerance','-Tolerance',...
                'Orientation','horizontal');
        p2 = subplot(2,1,2);
            h6 = plot(tout/3600,rad2deg(simulationsC(k).magnOmega),'-','LineWidth',2,'Color','#4DBEEE'); hold on; grid on;
            max = 3*0.1*ones(1,length(tout));
            h7 = plot(tout/3600,max_tolerance_norm,'k--','LineWidth',2);
            h8 = plot(simulationsC(k).tiempo_inicio/3600, 0.3, 'rx','MarkerSize',5);
            xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
            title('Simulations k constant: Satellite angular rate magnitude (deg/s)');
            legend([h6,h7,h8],'|\omega|','Tolerance','Settlement time','Orientation','horizontal');
            linkaxes([p1,p2],'x');
    else
        %%%Update angular rates
        h1xadata = get(h1,'XData'); h1ydata = get(h1,'YData');
        set(h1,'XData',[h1xadata,nan,tout'/3600],'YData',[h1ydata,nan,rad2deg(pqrout(:,1))']);
        %%%Update angular rates
        h2xadata = get(h2,'XData'); h2ydata = get(h2,'YData');
        set(h2,'XData',[h2xadata,nan,tout'/3600],'YData',[h2ydata,nan,rad2deg(pqrout(:,2))']);
        %%%Update angular rates
        h3xadata = get(h3,'XData'); h3ydata = get(h3,'YData');
        set(h3,'XData',[h3xadata,nan,tout'/3600],'YData',[h3ydata,nan,rad2deg(pqrout(:,3))']);
        %%%Update angular rates Magnitude
        h6xadata = get(h6,'XData'); h6ydata = get(h6,'YData');
        set(h6,'XData',[h6xadata,nan,tout'/3600],'YData',[h6ydata,nan,rad2deg(simulationsC(k).magnOmega)]);
        %%%Update settlement time
        h8xadata = get(h8,'XData'); h8ydata = get(h8,'YData');
        set(h8,'XData',[h8xadata,nan,simulationsC(k).tiempo_inicio/3600],'YData',[h8ydata,nan,0.3]);
    end
end
settlementTime = settlementTime/num;
plot([settlementTime/3600, settlementTime/3600], [0, 15^1/2], 'r--','LineWidth',2);
anotation_setTime = ['Average Settlement time: ', num2str(settlementTime/3600), ' h'];
x = [0.3,0.24];
y = [0.3,0.2];
annotation('textarrow',x,y,'String',anotation_setTime)

%% Plots
figure()
for k = 1:num
% Rename index
Torque_mag_v = simulationsC(k).T_control;   %Angular Rates
tout = simulationsC(k).tout;
m_mag_v=simulationsC(k).mu;
mag_currents=simulationsC(k).mag_currents;

% p1 = subplot(3,1,1);
%     h1 = plot(tout/3600,Torque_mag_v(1,:),'b-','LineWidth',2); hold on; grid on
%     h2 = plot(tout/3600,Torque_mag_v(2,:),'r-','LineWidth',2);
%     h3 = plot(tout/3600,Torque_mag_v(3,:),'g-','LineWidth',2);
%     legend([h1,h2,h3],'T_x','T_y','T_z','Location','eastoutside');
%     xlabel('Time (h)'); ylabel('Magnetic Torque (Nm)');
%     title('GWSAT magnetic Torque');
%     
% p2 = subplot(3,1,2);
%     h1 = plot(tout/3600,m_mag_v(1,:),'b-','LineWidth',2); hold on; grid on
%     h2 = plot(tout/3600,m_mag_v(2,:),'r-','LineWidth',2);
%     h3 = plot(tout/3600,m_mag_v(3,:),'g-','LineWidth',2);
%     legend([h1,h2,h3],'m_x','m_y','m_z','Location','eastoutside');
%     xlabel('Time (h)'); ylabel('Magnetic moment (Am^2)');
%     title('GWSAT magnetic moment');

% p3 = subplot(3,1,3);
    h1 = plot(tout/3600,mag_currents(1,:),'b-','LineWidth',2); hold on; grid on
    h2 = plot(tout/3600,mag_currents(2,:),'r-','LineWidth',2);
    h3 = plot(tout/3600,mag_currents(3,:),'g-','LineWidth',2);
    legend([h1,h2,h3],'I_x','I_y','I_z','Location','eastoutside');
    xlabel('Time (h)'); ylabel('Magnetorquer current (A)');
    title('Magnetorquers currents');
% linkaxes([p1,p2,p3],'x');
end
%% plot
totalPowerAvr=0;
figure()
for k = 1:num
    Px = simulationsC(k).powerX;
    Py = simulationsC(k).powerY;
    Pz = simulationsC(k).powerZ;
    totalPower=simulationsC(k).TotalPower;
    totalPowerAvr = totalPowerAvr+totalPower;
p1 = subplot(2,1,1);
    h1=plot(k,Px,'b.','LineWidth',2); hold on; grid on;
    h2=plot(k,Py,'r.','LineWidth',2);
    h3=plot(k,Pz,'g.','LineWidth',2); 
    legend([h1,h2,h3],'P_x','P_y','P_z','Location','eastoutside');
    xlabel('Simulation Number'); ylabel('Average Power(W)');
    title('Average Power per axis');
p2 = subplot(2,1,2);
    h4 = plot(k,totalPower,'b.','LineWidth',2); hold on; grid on;
    legend(h4,'P_{total}','Location','eastoutside');
    xlabel('Simulation Number'); ylabel('Average Total Power(W)');
    title('Average Power');
%linkaxes([p1,p2],'x');
end
totalPowerAvr = totalPowerAvr/num;
disp(totalPowerAvr);




%% Plots k variable
figure()
%%%Set toletances
max_tolerance      = 0.1*ones(1,length(simulationsV(1).tout));
max_tolerance_norm = 0.3*ones(1,length(simulationsV(1).tout));
%%%Calculate average settlement time 
settlementTime = 0;
for k = 1:num
    %%% Sum SettlementTime
    settlementTime = settlementTime + simulationsV(k).tiempo_inicio;
    %%% Rename index
    pqrout = simulationsV(k).x(:,5:7);      %% Angular Rates
    tout = simulationsV(k).tout;            %% time obtained
    
    if k == 1
        p1 = subplot(2,1,1);
            h1 = plot(tout/3600,rad2deg(pqrout(:,1)),'b-','LineWidth',2); hold on; grid on
            h2 = plot(tout/3600,rad2deg(pqrout(:,2)),'r-','LineWidth',2);
            h3 = plot(tout/3600,rad2deg(pqrout(:,3)),'g--','LineWidth',2);
            h4 = plot(tout/3600,max_tolerance,'k--','LineWidth',1);
            h5 = plot(tout/3600,-max_tolerance,'k--','LineWidth',1);
            xlabel('Time (h)'); ylabel('\omega (deg/s)');
            title('Simulations k constant: Satellite angular rate (deg/s)');
            legend([h1,h2,h3,h4,h5],'\omega_x','\omega_y','\omega_z','Tolerance','-Tolerance',...
                'Orientation','horizontal');
        p2 = subplot(2,1,2);
            h6 = plot(tout/3600,rad2deg(simulationsV(k).magnOmega),'-','LineWidth',2,'Color','#4DBEEE'); hold on; grid on;
            max = 3*0.1*ones(1,length(tout));
            h7 = plot(tout/3600,max_tolerance_norm,'k--','LineWidth',2);
            h8 = plot(simulationsV(k).tiempo_inicio/3600, 0.3, 'rx','MarkerSize',5);
            xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
            title('Simulations k constant: Satellite angular rate magnitude (deg/s)');
            legend([h6,h7,h8],'|\omega|','Tolerance','Settlement time','Orientation','horizontal');
            linkaxes([p1,p2],'x');
    else
        %%%Update angular rates
        h1xadata = get(h1,'XData'); h1ydata = get(h1,'YData');
        set(h1,'XData',[h1xadata,nan,tout'/3600],'YData',[h1ydata,nan,rad2deg(pqrout(:,1))']);
        %%%Update angular rates
        h2xadata = get(h2,'XData'); h2ydata = get(h2,'YData');
        set(h2,'XData',[h2xadata,nan,tout'/3600],'YData',[h2ydata,nan,rad2deg(pqrout(:,2))']);
        %%%Update angular rates
        h3xadata = get(h3,'XData'); h3ydata = get(h3,'YData');
        set(h3,'XData',[h3xadata,nan,tout'/3600],'YData',[h3ydata,nan,rad2deg(pqrout(:,3))']);
        %%%Update angular rates Magnitude
        h6xadata = get(h6,'XData'); h6ydata = get(h6,'YData');
        set(h6,'XData',[h6xadata,nan,tout'/3600],'YData',[h6ydata,nan,rad2deg(simulationsV(k).magnOmega)]);
        %%%Update settlement time
        h8xadata = get(h8,'XData'); h8ydata = get(h8,'YData');
        set(h8,'XData',[h8xadata,nan,simulationsV(k).tiempo_inicio/3600],'YData',[h8ydata,nan,0.3]);
    end
end
settlementTime = settlementTime/num;
plot([settlementTime/3600, settlementTime/3600], [0, 15^1/2], 'r--','LineWidth',2);
anotation_setTime = ['Average Settlement time: ', num2str(settlementTime/3600), ' h'];
x = [0.3,0.24];
y = [0.3,0.2];
annotation('textarrow',x,y,'String',anotation_setTime)

%% Plots
figure()
for k = 1:num
% Rename index
Torque_mag_v = simulationsV(k).T_control;   %Angular Rates
tout = simulationsV(k).tout;
m_mag_v=simulationsV(k).mu;
mag_currents=simulationsV(k).mag_currents;

% p1 = subplot(3,1,1);
%     h1 = plot(tout/3600,Torque_mag_v(1,:),'b-','LineWidth',2); hold on; grid on
%     h2 = plot(tout/3600,Torque_mag_v(2,:),'r-','LineWidth',2);
%     h3 = plot(tout/3600,Torque_mag_v(3,:),'g-','LineWidth',2);
%     legend([h1,h2,h3],'T_x','T_y','T_z','Location','eastoutside');
%     xlabel('Time (h)'); ylabel('Magnetic Torque (Nm)');
%     title('GWSAT magnetic Torque');
%     
% p2 = subplot(3,1,2);
%     h1 = plot(tout/3600,m_mag_v(1,:),'b-','LineWidth',2); hold on; grid on
%     h2 = plot(tout/3600,m_mag_v(2,:),'r-','LineWidth',2);
%     h3 = plot(tout/3600,m_mag_v(3,:),'g-','LineWidth',2);
%     legend([h1,h2,h3],'m_x','m_y','m_z','Location','eastoutside');
%     xlabel('Time (h)'); ylabel('Magnetic moment (Am^2)');
%     title('GWSAT magnetic moment');

% p3 = subplot(3,1,3);
    h1 = plot(tout/3600,mag_currents(1,:),'b-','LineWidth',2); hold on; grid on
    h2 = plot(tout/3600,mag_currents(2,:),'r-','LineWidth',2);
    h3 = plot(tout/3600,mag_currents(3,:),'g-','LineWidth',2);
    legend([h1,h2,h3],'I_x','I_y','I_z','Location','eastoutside');
    xlabel('Time (h)'); ylabel('Magnetorquer current (A)');
    title('Magnetorquers currents');
% linkaxes([p1,p2,p3],'x');
end

%% plot comparisson between sim1 and sim2
%%%Containers to calculate the power average
totalPowerAvr_const = 0;
totalPowerAvr_var   = 0;
%%%Initialize Containers
Px_const = zeros(1,num);Py_const = zeros(1,num); Pz_const = zeros(1,num);
Px_var = zeros(1,num);  Py_var = zeros(1,num); Pz_var = zeros(1,num);
totalPower_const = zeros(1,num); totalPower_var = zeros(1,num);

%%%Recover values
for k = 1:num
    Px_const(k) = simulationsC(k).powerX;
    Py_const(k) = simulationsC(k).powerY;
    Pz_const(k) = simulationsC(k).powerZ;
    
    Px_var(k)   = simulationsV(k).powerX;
    Py_var(k)   = simulationsV(k).powerY;
    Pz_var(k)   = simulationsV(k).powerZ;
    
    totalPower_const(k) = simulationsC(k).TotalPower;
    totalPower_var(k)   = simulationsV(k).TotalPower;
    
    totalPowerAvr_const = totalPowerAvr_const + totalPower_const(k);
    totalPowerAvr_var   = totalPowerAvr_var   + totalPower_var(k);
    
end
totalPowerAvr_const = totalPowerAvr_const/num;
totalPowerAvr_var   = totalPowerAvr_var  /num;

figure
   % p1 = subplot(2,1,1);
        h1 = plot(1:k,totalPower_const,'b-','LineWidth',2); hold on; grid on;
        h2 = plot(1:k,totalPower_var,'r-','LineWidth',2);
        %h3 = plot(1:k,totalPowerAvr_const*(ones(1,k)),'b--','LineWidth',2);
        %h4 = plot(1:k,totalPowerAvr_var*(ones(1,k)),'r--','LineWidth',2);
        %legend([h1,h2,h3,h4],'P_{kconst}','P_{kvar}','mean','mean','Location','eastoutside');
        legend([h1,h2],['P_{kconst}, mean: ' num2str(totalPowerAvr_const) 'W'],...
                        ['P_{kvar} mean: ' num2str(totalPowerAvr_var) 'W'],'Orientation','horizontal');
        xlabel('Simulation Number'); ylabel('Average Power(W)');
        title('Average Power');






%% Functions
function sim = detumblingSimulation(sim, kfunConstant, sat, disturbance, sensors, settings,actuators, X0)

%% Solve by variable step Runge Kutta of 4 order and use SPG4 propagator
%opts = odeset('InitialStep',1e-4,'RelTol',1e-6);
settings.hWaitbar = waitbar(0, 'Progress: 0%','Name', 'GWSAT-1 Simulation Progress');
%[tout,x] = ode45(@(t,x) satelliteSPG4sim(t, x, sat, disturbance, sensors, settings),...
%    [0 settings.tfinal],settings.X0,opts);

[tout,x] = ode45(@(t,x) satelliteSPG4sim(t, x, kfunConstant, sat, disturbance, sensors, settings),...
    [0 settings.tfinal],X0);

%%Asign variables
sim.tout = tout;
sim.x = x;

%% Recover values
%%% Create Data containers
T_control = zeros(3,length(tout));
mu        = zeros(3,length(tout));
mag_currents = zeros(3,length(tout));
magnOmega = zeros(1,length(tout));

settings.hWaitbar = waitbar(0, 'Progress: 0%','Name', 'GWSAT-1 Simulation Recover');
t = tout; 

for i = 1:length(x)
    %%Global variables
%global tant;

%% Disturbances
%LS2125204: Brayan Espinoza
%%%Calculate mag Omega
magnOmega(:,i) = norm(x(i,5:7)); 

%%% Update time vector
datetime_t = settings.startTime  + seconds(tout(i));

%%% Get state in geodetic and inertial reference frames.
[positionGCRF, ~] = states(sat.satSGP4, datetime_t, "CoordinateFrame", "inertial");
r_ecef = eci2ecef(datevec(datetime_t),positionGCRF);
llaSPG4 = ecef2lla(r_ecef', 'WGS84');
%%% Get magnetic field
[B_ref, ~, ~, dip, ~] = wrldmagm(llaSPG4(3), llaSPG4(1),llaSPG4(2),...
         settings.B_earth_field.decimal_year, settings.B_earth_field.model_epoch);  
k=kfunConstant(dip);
%%% Turn reference to BodyFrame
B_body = quatRotation(quatconj(x(i,1:4)), B_ref*1E-9);
%%% Apply magnetometer model
%%%mag_bm = mag_model(B_body,sensors.mag.desvEst,sensors.mag.res);

%%% Satellite model
%k = 2*(2*pi/5.677016087140827e+03)*(1+sin(deg2rad(dip)))*min(0.0789)*8e9;
[T_control(:,i),mu(:,i)] = detumblingControl(x(i,:)',k,B_body);
 %%% Obtain currents based on model
 mag_currents(:,i) = mu(:,i)/(actuators.magnetorquer.n*actuators.magnetorquer.A);
 
%%% Progress bar configuration
%%%% Calculate the progress percentage
progress = t(i) / settings.tfinal;

if isa(settings.hWaitbar,'handle') && isvalid(settings.hWaitbar)
    if t(i) == settings.tfinal
        % Close the progress bar when the simulation is complete
        close(settings.hWaitbar);
    else
        % Update the progress bar
        waitbar(progress, settings.hWaitbar, sprintf('Progress: %.1f%%', progress*100));
    end
end
end
sim.T_control=T_control;
sim.mu=mu;
sim.mag_currents=mag_currents;
sim.magnOmega=magnOmega;

%% Calculate parameters
[tiempo_inicio, index] = obtainStableTime(tout, magnOmega, 0, deg2rad(0.3));
powerX=potenciaMedia(tout(1:index)',mag_currents(1,1:index)); 
powerY=potenciaMedia(tout(1:index)',mag_currents(2,1:index));
powerZ=potenciaMedia(tout(1:index)',mag_currents(3,1:index));

TotalPower = powerX+powerY+powerZ;
sim.tiempo_inicio=tiempo_inicio;
sim.powerX = powerX;
sim.powerY = powerY;
sim.powerZ = powerZ;
sim.TotalPower = TotalPower;

end

%% functions
function [x_dot] = mySatelliteSPG4(t, x, sat, Td)
%%% Inputs
% t: current time
% x: satellite state
% sat: Satellite parameters
% earth: Earth parameters
% Td: External Disturbances

%%% Read state input
q   =  x(1:4);
w   =  x(5:7);

%% Read Satellite physical parameters
%m  = sat.mass;
Is = sat.Is;

%%%Rotational Dynamics
%%Book: Fundamentals of Spacecraft Attitude Determination and Control
%%Author: F. Landis Markley & John L. Crassidis
%%Quaternions Kinematics..Ecuación (2.88)
Xi=[-q(2),-q(3),-q(4);
    q(1),-q(4),q(3);
    q(4),q(1),-q(2);
    -q(3),q(2),q(1)];
q_dot=1/2*Xi*w;
%%%Kynematics Equation(3.21)
H = Is*w;                         %Satellite Angular momentun
w_dot = Is\(Td - cross(w,H));     %Dynamics Equatión (3.147)
%%%x_dot Vector
x_dot=[q_dot;w_dot];
end

function [x_dot] = satelliteSPG4sim(t, x, kfun, sat, dist, sensors, settings)
%%Global variables
%global tant;

%% Disturbances
%LS2125204: Brayan Espinoza

%%% Update time vector
datetime_t = settings.startTime  + seconds(t);

%%% Use Simplified disturbances
%T_disturbances = dist(t);

%%% Get state in geodetic and inertial reference frames.
[positionGCRF, ~] = states(sat.satSGP4, datetime_t, "CoordinateFrame", "inertial");
r_ecef = eci2ecef(datevec(datetime_t),positionGCRF);
llaSPG4 = ecef2lla(r_ecef', 'WGS84');
%%% Get magnetic field
[B_ref, ~, ~, dip, ~] = wrldmagm(llaSPG4(3), llaSPG4(1),llaSPG4(2),...
         settings.B_earth_field.decimal_year, settings.B_earth_field.model_epoch);  

%%% Turn reference to BodyFrame
B_body = quatRotation(quatconj(x(1:4)'), B_ref*1E-9);
%%% Apply magnetometer model
%%%mag_bm = mag_model(B_body,sensors.mag.desvEst,sensors.mag.res);

%%% Satellite model
%
k=kfun(dip);
[T_control,~] = detumblingControl(x,k,B_body);
%TotalTorque = T_disturbances+T_control';
%Sat function;
x_dot = mySatelliteSPG4(t, x, sat, T_control');

%%% Progress bar configuration
%%%% Calculate the progress percentage
progress = t / settings.tfinal;

if isa(settings.hWaitbar,'handle') && isvalid(settings.hWaitbar)
    if t == settings.tfinal
        % Close the progress bar when the simulation is complete
        close(settings.hWaitbar);
    else
        % Update the progress bar
        waitbar(progress, settings.hWaitbar, sprintf('Progress: %.1f%%', progress*100));
    end
end
end

function k=variableK(dip)
    k = 2*(2*pi/5.677016087140827e+03)*(1+sin(deg2rad(dip)))*min(0.0789)*8e9;
end

function k=constantK(dip)
    k = 0.25E6;
end
function rotX = quatRotation(q,x)
    qx = [0,x(1),x(2),x(3)];
    [a,b]=size(q);
    if a == 4 && b == 1
        q = q';
    end
    qrotX = quatmultiply(quatmultiply(q, qx), quatconj(q));
    %qrotX = QuaternionsMultiplication(QuaternionsMultiplication(q,qx),QuaternionsConjugate(q));
    rotX = qrotX(2:4);
end

function [Tc,muB] = detumblingControl(state,k,B_body)
 %% Earth magnetic field expresend in the body Frame
omega = state(5:7); 
muB = k*cross(omega,B_body);
Tc = cross(muB,B_body);
end

function [tiempo_inicio, index] = obtainStableTime(time, signal, lowerLimit, upperLimit)
    % Supongamos que 'time' es el vector de tiempo y 'signal' es la señal
    % 'lower' y 'upper' son los límites de la franja

    % Encuentra los índices donde la señal está dentro de la franja
    indices_en_franja = find(signal >= lowerLimit & signal <= upperLimit);

    % Verifica si la señal entra en la franja en algún momento
    if ~isempty(indices_en_franja)
        % Encuentra el tiempo correspondiente a los primeros y últimos índices en la franja
        index = indices_en_franja(1);
        tiempo_inicio = time(index);
        %tiempo_final = time(indices_en_franja(end));

        %% Calcula la duración en la franja
        %duracion_en_franja = tiempo_final - tiempo_inicio;

        %disp(['La señal está dentro de la franja desde ', num2str(tiempo_inicio), ' hasta ', num2str(tiempo_final)]);
        %disp(['La duración en la franja es ', num2str(duracion_en_franja), ' unidades de tiempo.']);
    else
        disp('La señal no está dentro de la franja en ningún momento.');
    end
end

function potencia_media_rectangular = potenciaMedia(t,intensity)
% Supongamos que tienes una señal de corriente en el vector i y el vector de tiempo t

% Calcula la potencia media utilizando la regla del punto medio
dt = diff(t);
%t_midpoints = t(1:end-1) + dt/2;
potencia_media_rectangular = sum(intensity(1:end-1).^2 .* dt) / (t(end) - t(1));

%disp(['Potencia Media (Método Rectangular): ', num2str(potencia_media_rectangular), ' vatios']);
end 