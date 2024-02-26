close all;  %%kvarkconstcomplete.mat
%% Satellite parameters
sat.Is = [0.103682,-0.001094, 0.015519;...
         -0.001094, 0.104158,-2.57e-6;...
          0.015519,  2.57e-6, 0.088588];
sat.mass = 40; %Kg

%% Actuators parameters
actuators.magnetorquer.power            = 1536E-3;       %W
actuators.magnetorquer.MaxPower         = 13200-3;       %W
actuators.magnetorquer.voltage          = 12;            %V
actuators.magnetorquer.dimensions       = [250,67,35];   %mm (1E-3 m)
actuators.magnetorquer.nominalDipole    = 15;            %Am^2
actuators.magnetorquer.maxNominalDipole = 30;            %Am^2
actuators.magnetorquer.A  = actuators.magnetorquer.dimensions(2)*actuators.magnetorquer.dimensions(3)*1E-6;   %Transversal Area (m^2)
actuators.magnetorquer.n  = actuators.magnetorquer.nominalDipole*actuators.magnetorquer.voltage/(actuators.magnetorquer.A*actuators.magnetorquer.power);           %Number of turns
actuators.magnetorquer.maxCurrent = actuators.magnetorquer.maxNominalDipole/actuators.magnetorquer.voltage;  %Max current (mA)

%% Sensor parameters
sensors.mag.desvEst = 50.00e-9;  %T
sensors.mag.res     = 31.25e-9;  %T
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

%% Function Gains
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

%% Start Random Simulations
num = 30; %%%Number of simulations for Monte Carlo analisys
%%%Simulations for B dot Controller 1
simulation1 = struct();
simulations1 = struct();
%%%Simulations for B-dot controller 2
simulation2 = struct();
simulations2 = struct();
X0 = zeros(7,num);

%% Generate random initial conditions
for k = 1:num
    X0(:,k) =[initial.atitude.q0123_0; deg2rad((rand(1) * 10)-5); ...
              deg2rad((rand(1) * 10)-5); deg2rad((rand(1) * 10)-5)]; 
end

%% Simulation for B-dot controller 1
for k = 1:num
    simulation1 = detumblingSimulation(simulation1,kfunConstant, ...
                       sat, disturbance, sensors, settings, actuators, X0(:,k));
    simulation1.X0 = X0(:,k);
    if k == 1
        simulations1 = simulation1;
        figure()
    else
        simulations1 = [simulations1;simulation1];
    end
    %% Plots
    % Rename index
    %q0123out = simulation1.x(:,1:4);              %%Satellite quaternions
    %ptpout = Quaternions2EulerAngles(q0123out);   %%Satellite euler angles
    pqrout = simulation1.x(:,5:7);                 %%Angular Rates
    tout = simulation1.tout;                       %%Simulation Time
    %Start simulation
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
        h1 = plot(tout/3600,rad2deg(simulation1.magnOmega),'LineWidth',2); hold on; grid on;
        max = 3*0.1*ones(1,length(tout));
        h2 = plot(tout/3600,max,'k--','LineWidth',2);
        xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
        title('Satellite angular rate magnitude (deg/s)');
        %xlim([0.0 5.0]);
        legend([h1,h2],'|\omega|','Tolerance','Location','eastoutside');
        linkaxes([p1,p2],'x')
end

%% Simulation for B-dot Controller 2
max_tolerance      = 0.1*ones(1,length(tout));
max_tolerance_norm = 0.3*ones(1,length(tout));
for k = 1:num
    simulation2 = detumblingSimulation(simulation2,kfunVariable, ...
                       sat, disturbance, sensors, settings, actuators, X0(:,k));
    simulation2.X0 = X0(:,k);
    if k == 1
        simulations2 = simulation2;
        
        figure()
        %%% Rename index
        pqrout = simulation2.x(:,5:7);      %%Angular Rates
        tout = simulation2.tout;            %%Vector time
        max_tolerance      = 0.1*ones(1,length(tout));
        max_tolerance_norm = 0.3*ones(1,length(tout));
        p1 = subplot(2,1,1);
            h1 = plot(tout/3600,rad2deg(pqrout(:,1)),'b-','LineWidth',2); hold on; grid on
            h2 = plot(tout/3600,rad2deg(pqrout(:,2)),'r-','LineWidth',2);
            h3 = plot(tout/3600,rad2deg(pqrout(:,3)),'g-','LineWidth',2);
            h4 = plot(tout/3600,max_tolerance,'k--','LineWidth',1);
            h5 = plot(tout/3600,-max_tolerance,'k--','LineWidth',1);
            xlabel('Time (h)'); ylabel('\omega (deg/s)');
            legend([h1,h2,h3,h4,h5],'\omega_x','\omega_y','\omega_z','Tolerance','Tolerance','Location','eastoutside');
            title('Satellite angular rate (deg/s)');
            %xlim([0.0 5.0]);

        p2 = subplot(2,1,2);
            h6 = plot(tout/3600,rad2deg(simulation2.magnOmega),'LineWidth',2); hold on; grid on;
            h7 = plot(tout/3600,max_tolerance_norm,'k--','LineWidth',2);
            xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
            title('Satellite angular rate magnitude (deg/s)');
            legend([h6,h7],'|\omega|','Tolerance','Location','eastoutside');
            linkaxes([p1,p2],'x');
    else
        simulations2 = [simulations2;simulation2];
    end
end

%% Plots initial conditions
figure()
plot(rad2deg(X0(5:7,:))','o'); hold on; grid on;
xlabel('Number of simulation'); ylabel('\omega_0 (deg/s)');
legend('\omega_{x0}', '\omega_{y0}', '\omega_{x0}');

%% Plot disturbances
figure()
plot(simulations1(1).tout/3600,simulations1(1).T_dist*1E3,'.');grid on;
xlabel('Time (h)'); ylabel('Torque (mNm)');
title('Disturbance torque in body Frame'); 
legend('T_x','T_y','T_z','Orientation','Horizontal');

%% Plot magnetic field
figure()
plot(simulations1(1).tout/3600,simulations1(1).mag_bm*1E9,'.');grid on;
xlabel('Time (h)'); ylabel('Magnetic field (nT)');
title('Magnetic field expressed in body frame in body Frame'); 
legend('Bm_x','Bm_y','Bm_z','Orientation','Horizontal');

%% Plots simulation 1
figure()
%%%Set toletances
max_tolerance      = 0.1*ones(1,length(simulations1(1).tout));
max_tolerance_norm = 0.3*ones(1,length(simulations1(1).tout));
%%%Calculate average settlement time 
settlementTime = 0;
for k = 1:num
    %%% Sum SettlementTime
    settlementTime = settlementTime + simulations1(k).tiempo_inicio;
    %%% Rename index
    pqrout = simulations1(k).x(:,5:7);      %% Angular Rates
    tout = simulations1(k).tout;            %% time obtained
    
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
            h6 = plot(tout/3600,rad2deg(simulations1(k).magnOmega),'-','LineWidth',2,'Color','#4DBEEE'); hold on; grid on;
            max = 3*0.1*ones(1,length(tout));
            h7 = plot(tout/3600,max_tolerance_norm,'k--','LineWidth',2);
            h8 = plot(simulations1(k).tiempo_inicio/3600, 0.3, 'rx','MarkerSize',5);
            xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
            title('Simulations k constant: Satellite angular rate magnitude (deg/s)');
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
        set(h6,'XData',[h6xadata,nan,tout'/3600],'YData',[h6ydata,nan,rad2deg(simulations1(k).magnOmega)]);
        %%%Update settlement time
        h8xadata = get(h8,'XData'); h8ydata = get(h8,'YData');
        set(h8,'XData',[h8xadata,nan,simulations1(k).tiempo_inicio/3600],'YData',[h8ydata,nan,0.3]);
    end
end
settlementTime = settlementTime/num;
h9 = plot([settlementTime/3600, settlementTime/3600], [0, 15^1/2], 'r--','LineWidth',2);
 legend([h6,h7,h8,h9],'|\omega|','Tolerance','Settlement time','Mean','Orientation','horizontal');
anotation_setTime = ['Average Settlement time: ', num2str(settlementTime/3600), ' h'];
x = [0.3,0.24];
y = [0.3,0.2];
annotation('textarrow',x,y,'String',anotation_setTime)

%% Plot simulation 1 remain torques
figure()
%%%Calculate average settlement time 
settlementTime = 0;
for k = 1:num
    %%% Sum SettlementTime
    settlementTime = settlementTime + simulations1(k).tiempo_inicio;
    %%% Rename index
    pqrout = simulations1(k).x(:,5:7);      %% Angular Rates
    tout = simulations1(k).tout;            %% time obtained
    
    if k == 1
        h6 = plot(tout/3600,rad2deg(simulations1(k).magnOmega),'-','LineWidth',2,'Color','#4DBEEE'); hold on; grid on;
        max = 3*0.1*ones(1,length(tout));
        h7 = plot(tout/3600,max_tolerance_norm,'k--','LineWidth',2);
        h8 = plot(simulations1(k).tiempo_inicio/3600, 0.3, 'rx','MarkerSize',5);
        xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
        title('Simulations k constant: Satellite angular rate magnitude (deg/s)');
    else
        %%%Update angular rates Magnitude
        h6xadata = get(h6,'XData'); h6ydata = get(h6,'YData');
        set(h6,'XData',[h6xadata,nan,tout'/3600],'YData',[h6ydata,nan,rad2deg(simulations1(k).magnOmega)]);
        %%%Update settlement time
        h8xadata = get(h8,'XData'); h8ydata = get(h8,'YData');
        set(h8,'XData',[h8xadata,nan,simulations1(k).tiempo_inicio/3600],'YData',[h8ydata,nan,0.3]);
    end
end
settlementTime = settlementTime/num;
h9 = plot([settlementTime/3600, settlementTime/3600], [0, 15^1/2], 'r--','LineWidth',2);
 legend([h6,h7,h8,h9],'|\omega|','Tolerance','Settlement time','Mean','Orientation','horizontal');
anotation_setTime = ['Average Settlement time: ', num2str(settlementTime/3600), ' h'];
x = [0.3,0.24];
y = [0.3,0.2];
annotation('textarrow',x,y,'String',anotation_setTime)
ylim([0,0.4]);

%% Plots simulation 1 current and magnetic moment
figure()
for k = 1:num
    % Rename index
    Torque_mag_v = simulations1(k).T_control;   %Torque delivered by magnetorquers
    tout = simulations1(k).tout;                %Recover simulation time
    m_mag_v=simulations1(k).mu;                 %Recover magnetic momentum
    mag_currents=simulations1(k).mag_currents;  %Recover magnetic currents
    
    if(k==1)
    p1 = subplot(3,1,1);
        h1 = plot(tout/3600,Torque_mag_v(1,:)*1E3,'b.','LineWidth',2); hold on; grid on
        h2 = plot(tout/3600,Torque_mag_v(2,:)*1E3,'r.','LineWidth',2);
        h3 = plot(tout/3600,Torque_mag_v(3,:)*1E3,'g.','LineWidth',2);
        legend([h1,h2,h3],'T_x','T_y','T_z','Orientation','horizontal');
        xlabel('Time (h)'); ylabel('T_m (mNm)');
        title('Simulations k constant: GWSAT magnetic Torque');
        
    p2 = subplot(3,1,2);
        h4 = plot(tout/3600,m_mag_v(1,:),'b.','LineWidth',2); hold on; grid on
        h5 = plot(tout/3600,m_mag_v(2,:),'r.','LineWidth',2);
        h6 = plot(tout/3600,m_mag_v(3,:),'g.','LineWidth',2);
        legend([h4,h5,h6],'m_x','m_y','m_z','Orientation','horizontal');
        xlabel('Time (h)'); ylabel('m_{mag} (Am^2)');
        title('Simulations k constant: GWSAT magnetic moment');

    p3 = subplot(3,1,3);
        h7 = plot(tout/3600,mag_currents(1,:)*1E3,'b.','LineWidth',2); hold on; grid on
        h8 = plot(tout/3600,mag_currents(2,:)*1E3,'r.','LineWidth',2);
        h9 = plot(tout/3600,mag_currents(3,:)*1E3,'g.','LineWidth',2);
        legend([h7,h8,h9],'I_x','I_y','I_z','Orientation','horizontal');
        xlabel('Time (h)'); ylabel('I_m (mA)');
        title('Simulations k constant: Magnetorquers currents');
    linkaxes([p1,p2,p3],'x');
    else
        %%%Update torquesx 
        h1xadata = get(h1,'XData'); h1ydata = get(h1,'YData');
        set(h1,'XData',[h1xadata,nan,tout'/3600],'YData',[h1ydata,nan,Torque_mag_v(1,:)*1E3]);
        %%%Update torquey
        h2xadata = get(h2,'XData'); h2ydata = get(h2,'YData');
        set(h2,'XData',[h2xadata,nan,tout'/3600],'YData',[h2ydata,nan,Torque_mag_v(2,:)*1E3]);
        %%%Update torquez
        h3xadata = get(h3,'XData'); h3ydata = get(h3,'YData');
        set(h3,'XData',[h3xadata,nan,tout'/3600],'YData',[h3ydata,nan,Torque_mag_v(3,:)*1E3]);
        
        %%%Update tmagmomentX  
        h4xadata = get(h4,'XData'); h4ydata = get(h4,'YData');
        set(h4,'XData',[h4xadata,nan,tout'/3600],'YData',[h4ydata,nan,m_mag_v(1,:)]);
        %%%Update tmagmomentY 
        h5xadata = get(h5,'XData'); h5ydata = get(h5,'YData');
        set(h5,'XData',[h5xadata,nan,tout'/3600],'YData',[h5ydata,nan,m_mag_v(2,:)]);
        %%%Update tmagmomentZ
        h6xadata = get(h6,'XData'); h6ydata = get(h6,'YData');
        set(h6,'XData',[h6xadata,nan,tout'/3600],'YData',[h6ydata,nan,m_mag_v(3,:)]);
        
        %%%Update currentx 
        h7xadata = get(h7,'XData'); h7ydata = get(h7,'YData');
        set(h7,'XData',[h7xadata,nan,tout'/3600],'YData',[h7ydata,nan,mag_currents(1,:)*1E3]);
        %%%Update currenty
        h8xadata = get(h8,'XData'); h8ydata = get(h8,'YData');
        set(h8,'XData',[h8xadata,nan,tout'/3600],'YData',[h8ydata,nan,mag_currents(2,:)*1E3]);
        %%%Update currentz
        h9xadata = get(h9,'XData'); h9ydata = get(h9,'YData');
        set(h9,'XData',[h9xadata,nan,tout'/3600],'YData',[h9ydata,nan,mag_currents(3,:)*1E3]);
    end
end

%% Plots simulation 2
figure()
%%%Set toletances
max_tolerance      = 0.1*ones(1,length(simulations2(1).tout));
max_tolerance_norm = 0.3*ones(1,length(simulations2(1).tout));
%%%Calculate average settlement time 
settlementTime = 0;
for k = 1:num
    %%% Sum SettlementTime
    settlementTime = settlementTime + simulations2(k).tiempo_inicio;
    %%% Rename index
    pqrout = simulations2(k).x(:,5:7);      %% Angular Rates
    tout = simulations2(k).tout;            %% time obtained
    
    if k == 1
        p1 = subplot(2,1,1);
            h1 = plot(tout/3600,rad2deg(pqrout(:,1)),'b-','LineWidth',2); hold on; grid on
            h2 = plot(tout/3600,rad2deg(pqrout(:,2)),'r-','LineWidth',2);
            h3 = plot(tout/3600,rad2deg(pqrout(:,3)),'g--','LineWidth',2);
            h4 = plot(tout/3600,max_tolerance,'k--','LineWidth',1);
            h5 = plot(tout/3600,-max_tolerance,'k--','LineWidth',1);
            xlabel('Time (h)'); ylabel('\omega (deg/s)');
            title('Simulations k variable: Satellite angular rate (deg/s)');
            legend([h1,h2,h3,h4,h5],'\omega_x','\omega_y','\omega_z','Tolerance','-Tolerance',...
                'Orientation','horizontal');
        p2 = subplot(2,1,2);
            h6 = plot(tout/3600,rad2deg(simulations2(k).magnOmega),'-','LineWidth',2,'Color','#4DBEEE'); hold on; grid on;
            max = 3*0.1*ones(1,length(tout));
            h7 = plot(tout/3600,max_tolerance_norm,'k--','LineWidth',2);
            h8 = plot(simulations2(k).tiempo_inicio/3600, 0.3, 'rx','MarkerSize',5);
            xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
            title('Simulations k variable: Satellite angular rate magnitude (deg/s)');
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
        set(h6,'XData',[h6xadata,nan,tout'/3600],'YData',[h6ydata,nan,rad2deg(simulations2(k).magnOmega)]);
        %%%Update settlement time
        h8xadata = get(h8,'XData'); h8ydata = get(h8,'YData');
        set(h8,'XData',[h8xadata,nan,simulations2(k).tiempo_inicio/3600],'YData',[h8ydata,nan,0.3]);
    end
end
settlementTime = settlementTime/num;
plot([settlementTime/3600, settlementTime/3600], [0, 15^1/2], 'r--','LineWidth',2);
anotation_setTime = ['Average Settlement time: ', num2str(settlementTime/3600), ' h'];
x = [0.3,0.24];
y = [0.3,0.2];
annotation('textarrow',x,y,'String',anotation_setTime)

%% Plot simulation 2 remain torques
figure()
%%%Calculate average settlement time 
settlementTime = 0;
max_tolerance_norm = 0.3*ones(1,length(simulations2(1).tout));

for k = 1:num
    %%% Sum SettlementTime
    settlementTime = settlementTime + simulations2(k).tiempo_inicio;
    %%% Rename index
    pqrout = simulations2(k).x(:,5:7);      %% Angular Rates
    tout = simulations2(k).tout;            %% time obtained
    
    if k == 1
        h6 = plot(tout/3600,rad2deg(simulations2(k).magnOmega),'-','LineWidth',2,'Color','#4DBEEE'); hold on; grid on;
        max = 3*0.1*ones(1,length(tout));
        h7 = plot(tout/3600,max_tolerance_norm,'k--','LineWidth',2);
        h8 = plot(simulations2(k).tiempo_inicio/3600, 0.3, 'rx','MarkerSize',5);
        xlabel('Time (h)'); ylabel('|\omega| (deg/s)');
        title('Simulations k variable: Satellite angular rate magnitude (deg/s)');
    else
        %%%Update angular rates Magnitude
        h6xadata = get(h6,'XData'); h6ydata = get(h6,'YData');
        set(h6,'XData',[h6xadata,nan,tout'/3600],'YData',[h6ydata,nan,rad2deg(simulations2(k).magnOmega)]);
        %%%Update settlement time
        h8xadata = get(h8,'XData'); h8ydata = get(h8,'YData');
        set(h8,'XData',[h8xadata,nan,simulations2(k).tiempo_inicio/3600],'YData',[h8ydata,nan,0.3]);
    end
end
settlementTime = settlementTime/num;
h9 = plot([settlementTime/3600, settlementTime/3600], [0, 15^1/2], 'r--','LineWidth',2);
 legend([h6,h7,h8,h9],'|\omega|','Tolerance','Settlement time','Mean','Orientation','horizontal');
anotation_setTime = ['Average Settlement time: ', num2str(settlementTime/3600), ' h'];
x = [0.3,0.24];
y = [0.3,0.2];
annotation('textarrow',x,y,'String',anotation_setTime)
ylim([0,0.4]);

%% Plots simulation 2 power
figure()
for k = 1:num
    % Rename index
    Torque_mag_v = simulations2(k).T_control;   %Torque delivered by magnetorquers
    tout = simulations2(k).tout;                %Recover simulation time
    m_mag_v=simulations2(k).mu;                 %Recover magnetic momentum
    mag_currents=simulations2(k).mag_currents;  %Recover magnetic currents
    
    if(k==1)
    p1 = subplot(3,1,1);
        h1 = plot(tout/3600,Torque_mag_v(1,:)*1E3,'b.','LineWidth',2); hold on; grid on
        h2 = plot(tout/3600,Torque_mag_v(2,:)*1E3,'r.','LineWidth',2);
        h3 = plot(tout/3600,Torque_mag_v(3,:)*1E3,'g.','LineWidth',2);
        legend([h1,h2,h3],'T_x','T_y','T_z','Orientation','horizontal');
        xlabel('Time (h)'); ylabel('T_m (mNm)');
        title('Simulations k variable: GWSAT magnetic Torque');
        
    p2 = subplot(3,1,2);
        h4 = plot(tout/3600,m_mag_v(1,:),'b.','LineWidth',2); hold on; grid on
        h5 = plot(tout/3600,m_mag_v(2,:),'r.','LineWidth',2);
        h6 = plot(tout/3600,m_mag_v(3,:),'g.','LineWidth',2);
        legend([h4,h5,h6],'m_x','m_y','m_z','Orientation','horizontal');
        xlabel('Time (h)'); ylabel('m_{mag} (Am^2)');
        title('Simulations k variable: GWSAT magnetic moment');

    p3 = subplot(3,1,3);
        h7 = plot(tout/3600,mag_currents(1,:)*1E3,'b.','LineWidth',2); hold on; grid on
        h8 = plot(tout/3600,mag_currents(2,:)*1E3,'r.','LineWidth',2);
        h9 = plot(tout/3600,mag_currents(3,:)*1E3,'g.','LineWidth',2);
        legend([h7,h8,h9],'I_x','I_y','I_z','Orientation','horizontal');
        xlabel('Time (h)'); ylabel('I_m (mA)');
        title('Simulations k variable: Magnetorquers currents');
    linkaxes([p1,p2,p3],'x');
    else
        %%%Update torquesx 
        h1xadata = get(h1,'XData'); h1ydata = get(h1,'YData');
        set(h1,'XData',[h1xadata,nan,tout'/3600],'YData',[h1ydata,nan,Torque_mag_v(1,:)*1E3]);
        %%%Update torquey
        h2xadata = get(h2,'XData'); h2ydata = get(h2,'YData');
        set(h2,'XData',[h2xadata,nan,tout'/3600],'YData',[h2ydata,nan,Torque_mag_v(2,:)*1E3]);
        %%%Update torquez
        h3xadata = get(h3,'XData'); h3ydata = get(h3,'YData');
        set(h3,'XData',[h3xadata,nan,tout'/3600],'YData',[h3ydata,nan,Torque_mag_v(3,:)*1E3]);
        
        %%%Update tmagmomentX  
        h4xadata = get(h4,'XData'); h4ydata = get(h4,'YData');
        set(h4,'XData',[h4xadata,nan,tout'/3600],'YData',[h4ydata,nan,m_mag_v(1,:)]);
        %%%Update tmagmomentY 
        h5xadata = get(h5,'XData'); h5ydata = get(h5,'YData');
        set(h5,'XData',[h5xadata,nan,tout'/3600],'YData',[h5ydata,nan,m_mag_v(2,:)]);
        %%%Update tmagmomentZ
        h6xadata = get(h6,'XData'); h6ydata = get(h6,'YData');
        set(h6,'XData',[h6xadata,nan,tout'/3600],'YData',[h6ydata,nan,m_mag_v(3,:)]);
        
        %%%Update currentx 
        h7xadata = get(h7,'XData'); h7ydata = get(h7,'YData');
        set(h7,'XData',[h7xadata,nan,tout'/3600],'YData',[h7ydata,nan,mag_currents(1,:)*1E3]);
        %%%Update currenty
        h8xadata = get(h8,'XData'); h8ydata = get(h8,'YData');
        set(h8,'XData',[h8xadata,nan,tout'/3600],'YData',[h8ydata,nan,mag_currents(2,:)*1E3]);
        %%%Update currentz
        h9xadata = get(h9,'XData'); h9ydata = get(h9,'YData');
        set(h9,'XData',[h9xadata,nan,tout'/3600],'YData',[h9ydata,nan,mag_currents(3,:)*1E3]);
    end
    
end

%% Plots simulation 2 gain
figure()
for k = 1:num
    % Rename index
    k_gain = simulations2(k).k;   %Recover gain
    tout = simulations2(k).tout;  %Recover simulation time
    if(k==1)
        h1 = plot(tout/3600,k_gain,'b.','LineWidth',2); grid on
        xlabel('Time (h)'); ylabel('Gain');
        title('Simulations k variable: B-dot gain');
    else
        %%%Update torquesx 
        h1xadata = get(h1,'XData'); h1ydata = get(h1,'YData');
        set(h1,'XData',[h1xadata,nan,tout'/3600],'YData',[h1ydata,nan,k_gain]);
    end
    
end

%% plot comparisson between sim1 and sim2
for k = 1:num
    % Rename index
    tout_1 = simulations1(k).tout;                  %Recover simulation time
    mag_currents_1 = simulations1(k).mag_currents;  %Recover magnetic currents
    
    tout_2 = simulations2(k).tout;                  %Recover simulation time
    mag_currents_2 = simulations2(k).mag_currents;  %Recover magnetic currents
    
    if(k==1)
    p1 = subplot(2,1,1);
        h1 = plot(tout_1/3600,mag_currents_1(1,:)*1E3,'r-','LineWidth',2); grid on; hold on
        h2 = plot(tout_1/3600,mag_currents_1(2,:)*1E3,'g-','LineWidth',2);
        h3 = plot(tout_1/3600,mag_currents_1(3,:)*1E3,'b-','LineWidth',2);
        legend([h1,h2,h3],'I_x','I_y','I_z','Orientation','horizontal');
        xlabel('Time (h)'); ylabel('I (mA)');
        title('Simulations k constant: Magnetorquers currents');
        
    p2 = subplot(2,1,2);
        h4 = plot(tout_2/3600,mag_currents_2(1,:)*1E3,'r-','LineWidth',2); grid on; hold on
        h5 = plot(tout_2/3600,mag_currents_2(2,:)*1E3,'g-','LineWidth',2); 
        h6 = plot(tout_2/3600,mag_currents_2(3,:)*1E3,'b-','LineWidth',2); 
        legend([h4,h5,h6],'I_x','I_y','I_z','Orientation','horizontal');
        xlabel('Time (h)'); ylabel('I (mA)');
        title('Simulations k variable: Magnetorquers currents');
    linkaxes([p1,p2],'x');
    
    else
        %%%Update ix
        h1xadata = get(h1,'XData'); h1ydata = get(h1,'YData');
        set(h1,'XData',[h1xadata,nan,tout_1'/3600],'YData',[h1ydata,nan,mag_currents_1(1,:)*1E3]);
        %%%Update iy
        h2xadata = get(h2,'XData'); h2ydata = get(h2,'YData');
        set(h2,'XData',[h2xadata,nan,tout_1'/3600],'YData',[h2ydata,nan,mag_currents_1(2,:)*1E3]);
        %%%Update iz
        h3xadata = get(h3,'XData'); h3ydata = get(h3,'YData');
        set(h3,'XData',[h3xadata,nan,tout_1'/3600],'YData',[h3ydata,nan,mag_currents_1(3,:)*1E3]);
        
        %%%Update ix
        h4xadata = get(h4,'XData'); h4ydata = get(h4,'YData');
        set(h4,'XData',[h4xadata,nan,tout_2'/3600],'YData',[h4ydata,nan,mag_currents_2(1,:)*1E3]);
        %%%Update iy
        h5xadata = get(h5,'XData'); h5ydata = get(h5,'YData');
        set(h5,'XData',[h5xadata,nan,tout_2'/3600],'YData',[h5ydata,nan,mag_currents_2(2,:)*1E3]);
        %%%Update iz
        h6xadata = get(h6,'XData'); h6ydata = get(h6,'YData');
        set(h6,'XData',[h6xadata,nan,tout_2'/3600],'YData',[h6ydata,nan,mag_currents_2(3,:)*1E3]);
    end
    
end

%% Plot calculated power
%%%Containers to calculate the power average
totalPowerAvr_const = 0;
totalPowerAvr_var   = 0;
%%%Initialize Containers
Px_const = zeros(1,num);Py_const = zeros(1,num); Pz_const = zeros(1,num);
Px_var = zeros(1,num);  Py_var = zeros(1,num); Pz_var = zeros(1,num);
totalPower_const = zeros(1,num); totalPower_var = zeros(1,num);

%%%Recover values
for k = 1:num
    Px_const(k) = simulations1(k).powerX;
    Py_const(k) = simulations1(k).powerY;
    Pz_const(k) = simulations1(k).powerZ;
    
    Px_var(k)   = simulations2(k).powerX;
    Py_var(k)   = simulations2(k).powerY;
    Pz_var(k)   = simulations2(k).powerZ;
    
    totalPower_const(k) = simulations1(k).TotalPower;
    totalPower_var(k)   = simulations2(k).TotalPower;
    
    totalPowerAvr_const = totalPowerAvr_const + totalPower_const(k);
    totalPowerAvr_var   = totalPowerAvr_var   + totalPower_var(k);
    
end
totalPowerAvr_const = totalPowerAvr_const/num;
totalPowerAvr_var   = totalPowerAvr_var  /num;

figure
    h1 = plot(1:k,totalPower_const,'b-','LineWidth',2); hold on; grid on;
    h2 = plot(1:k,totalPower_var,'r-','LineWidth',2);
    legend([h1,h2],['P_{kconst}, mean: ' num2str(totalPowerAvr_const) 'W'],...
            ['P_{kvar} mean: ' num2str(totalPowerAvr_var) 'W'],'Orientation','horizontal');
    xlabel('Simulation Number'); ylabel('Average Power(W)');
    title('Average Power');



%% Functions
function sim = detumblingSimulation(sim, kfunConstant, sat, disturbance, sensors, settings,actuators, X0)
    %% Initiate simulations
    %%% Call Loading bar
    settings.hWaitbar = waitbar(0, 'Progress: 0%','Name', 'GWSAT-1 Simulation Progress');
    %%% Solve by variable step Runge Kutta of 4 order and use SPG4 propagator
    [sim.tout,sim.x] = ode45(@(t,x) satelliteSPG4sim(t, x, kfunConstant, sat, disturbance, sensors, settings),...
        [0 settings.tfinal],X0);
    %% Recover values
    %%% Call Loading bar
    settings.hWaitbar = waitbar(0, 'Progress: 0%','Name', 'GWSAT-1 Simulation Progress');
    %%% Create Data containers
    mu           = zeros(3,length(sim.tout));
    T_control    = zeros(3,length(sim.tout));
    magnOmega    = zeros(1,length(sim.tout));
    mag_currents = zeros(3,length(sim.tout));
    kv           = zeros(1,length(sim.tout));
    mag_bm       = zeros(3,length(sim.tout));
    b_body       = zeros(3,length(sim.tout));
    b_ref        = zeros(3,length(sim.tout));
    T_dist       = zeros(3,length(sim.tout));
    
    t = sim.tout; 

    for i = 1:length(sim.x)
        %%Global variables
        %global tant;

        %% Disturbances
        %LS2125204: Brayan Espinoza
        %%%Calculate mag Omega
        magnOmega(:,i) = norm(sim.x(i,5:7)); 
        T_dist(:,i) = disturbance(t(i));
        %%% Update time vector
        datetime_t = settings.startTime  + seconds(sim.tout(i));

        %%% Get state in geodetic and inertial reference frames.
        [positionGCRF, ~] = states(sat.satSGP4, datetime_t, "CoordinateFrame", "inertial");
        r_ecef = eci2ecef(datevec(datetime_t),positionGCRF);
        llaSPG4 = ecef2lla(r_ecef', 'WGS84');
        %%% Get magnetic field
        [b_ref(:,i), ~, ~, dip, ~] = wrldmagm(llaSPG4(3), llaSPG4(1),llaSPG4(2),...
                 settings.B_earth_field.decimal_year, settings.B_earth_field.model_epoch);  
        kv(i)=kfunConstant(dip);
        %%% Turn reference to BodyFrame
        b_body(:,i) = quatRotation(quatconj(sim.x(i,1:4)), b_ref(:,i)*1E-9);
        %%% Apply magnetometer model
        mag_bm(:,i) = mag_model(b_body(:,i),sensors.mag.desvEst,sensors.mag.res);
        
        %%% Satellite model
        [T_control(:,i),mu(:,i)] = detumblingControl(sim.x(i,:)',kv(i),mag_bm(:,i));
         %%% Obtain currents based on model
         mag_currents(:,i) = mu(:,i)/(actuators.magnetorquer.n*actuators.magnetorquer.A);

        %%% Progress bar configuration
        %%%% Calculate the progress percentage
        progress = t(i) / settings.tfinal;

        if isa(settings.hWaitbar,'handle') && isvalid(settings.hWaitbar)
            if t(i) >= settings.tfinal
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
    sim.k = kv;
    sim.mag_bm = mag_bm;
    sim.b_body = b_body;
    sim.b_ref = b_ref;        
    sim.T_dist = T_dist;
    
    %% Calculate parameters
    [tiempo_inicio, index] = obtainStableTime(sim.tout, magnOmega, 0, deg2rad(0.3));
    powerX=potenciaMedia(sim.tout(1:index)',mag_currents(1,1:index)); 
    powerY=potenciaMedia(sim.tout(1:index)',mag_currents(2,1:index));
    powerZ=potenciaMedia(sim.tout(1:index)',mag_currents(3,1:index));

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
%% Disturbances
%LS2125204: Brayan Espinoza

    %%% Update time vector
    datetime_t = settings.startTime  + seconds(t);

    %%% Use Simplified disturbances
    T_disturbances = dist(t);

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
    mag_bm = mag_model(B_body,sensors.mag.desvEst,sensors.mag.res);

    %%% Get gain value
    k=kfun(dip);
    %%% Apply detumbling Controll
    [T_control,~] = detumblingControl(x,k,mag_bm);
    %TotalTorque = T_control';
    TotalTorque = T_disturbances+T_control';
    
    %Sat function;
    x_dot = mySatelliteSPG4(t, x, sat, TotalTorque);

    %% Progress bar configuration
    %%%% Calculate the progress percentage
    progress = t / settings.tfinal;

    if isa(settings.hWaitbar,'handle') && isvalid(settings.hWaitbar)
        if t >= settings.tfinal-1
            % Close the progress bar when the simulation is complete
            close(settings.hWaitbar);
        else
            % Update the progress bar
            waitbar(progress, settings.hWaitbar, sprintf('Progress: %.1f%%', progress*100));
        end
    end
end

function k=variableK(dip)
    %k = 2*(2*pi/5.677016087140827e+03)*(1+sin(deg2rad(dip)))*min(0.0789)*8e9;
    k = 2*(2*pi/5.677016087140827e+03)*(1+sin(deg2rad(dip)))*min(0.0789);
end

function k=constantK(dip)
    %k = 0.25E6;
    %k = 3.5E-4;
    %k = 0.25E6;
    k = 2.5E-4;
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
    b=B_body/norm(B_body);
    muB = k/norm(B_body)*cross(omega,b);
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

function mag_bm = mag_model(mag_b,sig_tam,sig_res)
    %%-----------magnetometerc Model------------
    %%%Add random noise
    mag_b_noise=mag_b+sig_tam*randn(1,3);
    %%%Add resolution
    mag_bm=[round(mag_b_noise(1)/sig_res)*sig_res,...
            round(mag_b_noise(2)/sig_res)*sig_res,...
            round(mag_b_noise(3)/sig_res)*sig_res];   
end