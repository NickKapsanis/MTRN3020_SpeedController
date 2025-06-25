% function to define constants into workspace for use in simulink 
%define all the physical constants
t_s = 0.009; %sampleing time
R_a = 4.89; %Ohm
L = 0.00042; %H
J = 0.0000109; % kg/m^2
B = 0.0000464; %Nm/(rad/s)
k_t = 0.0348; %Nm/A
k_e = 0.0348; %V/(rad/s)

%find sim time to run for
t = (1600-1)*t_s;

%define a few conversion factors
RPM2Cs = 8192/60;
Cs2RPM = 1/RPM2Cs;
Cs2Rads = 8192/(2*pi);
PWM2V = 24/126;

% get the transfer functions
Gc_z = Lab2_TransferFunction(); % this will create the Gp_z tf in the workspace

%run simulations
% load gain, if lg = 0, no load applied to sim, if lg = 1, step load applied to sim per zid.
lg = 0; % no load simulation
%run simulations
simA = sim('NoLoadSim', 'StopTime', num2str(t)); %no load simulation

lg = 1;
simB = sim('NoLoadSim', 'StopTime', num2str(t)); % Stepped load simulation

% Extract simulation results for no-load (A)
simA_speed = simA.shaftSpeed; % speed in counts/s
simA_set = simA.setSpeed; % set speed in counts/s

% Extract simulation results for no-load (B)
simB_speed = simB.shaftSpeed; % speed in counts/s
simB_set = simB.setSpeed; % set speed in counts/s

% get experimental data for part A, no load, and B, stepped load
run('C:\Users\nickk\OneDrive\Desktop\MTRN3020\LAB2_SpeedController\LAB2\Data\A5254990.m');
%run('C:\Users\nickk\OneDrive\Desktop\MTRN3020\LAB2_SpeedController\Lab2\Data\A5254990.M');
ExpA = RUN;
run('C:\Users\nickk\OneDrive\Desktop\MTRN3020\LAB2_SpeedController\Lab2\Data\B5254990.m');
ExpB = RUN;
clear RUN;

% Extract experimental data for no-load (ExpA)
timeA = ExpA(:, 1) / 1000; % Convert time from ms to seconds
expA_speed = ExpA(:, 3); % Shaft speed in counts/s
expA_set = ExpA(:, 4);   % Set shaft speed in counts/s

% Extract experimental data for stepped load (ExpB)
timeB = ExpB(:, 1) / 1000; % Convert time from ms to seconds
expB_speed = ExpB(:, 3); % Shaft speed in counts/s
expB_set = ExpB(:, 4);   % Set shaft speed in counts/s

% Plot the no-load data (ExpA)
figure;
plot(timeA, expA_set*Cs2RPM, 'r', 'LineWidth', 1.5); hold on; % Experimental set_speed in red
plot(timeA, expA_speed*Cs2RPM, 'b', 'LineWidth', 1.5); % Experimental shaft speed in blue
plot(timeA, simA_speed*Cs2RPM, 'g--', 'LineWidth', 1.5); % Simulated shaft speed in green dashed line
plot(timeA, simA_set*Cs2RPM, 'm--', 'LineWidth', 1.5); % Simulated set speed in magenta dashed line
xlabel('Time (s)');
ylabel('Speed (RPM)');
title('No-Load Experimental-Sim Comparison');
legend('Set Speed Experimental', 'Shaft Speed Experimental', 'Sim Shaft Speed', 'Sim Set-Speed');
grid on;

% Plot the load data (ExpB)
figure;
plot(timeB, expB_set*Cs2RPM, 'r', 'LineWidth', 1.5); hold on; % Experimental set_speed in red
plot(timeB, expB_speed*Cs2RPM, 'b', 'LineWidth', 1.5); % Experimental shaft speed in blue
plot(timeB, simB_speed*Cs2RPM, 'g--', 'LineWidth', 1.5); % Simulated shaft speed in green dashed line
plot(timeB, simB_set*Cs2RPM, 'm--', 'LineWidth', 1.5); % Simulated set speed in magenta dashed line
xlabel('Time (s)');
ylabel('Speed (RPM)');
title('Load Experimental-Sim Comparison');
legend('Set Speed Experimental', 'Shaft Speed Experimental', 'Sim Shaft Speed', 'Sim Set-Speed');
grid on;

% now analysis to prove tau_d = 36ms is correct
% Define the step time and time range
step_time = 2.07; % Time when the step occurs
time_limit = 0.3; % Stop plotting 0.3 seconds after the step - arbitrary but makes a pretty plot
tau = 0.036; % Desired time constant (36 ms)

% indices for to the step time and time limit
start_index = find(timeA >= step_time, 1); % Start when the 1000->2000 RPM step happens
end_index = find(timeA >= step_time + time_limit, 1);

% get data
timeA_step = timeA(start_index:end_index) - step_time; % rel time so that plotting shows settling time better
expA_speed_step = expA_speed(start_index:end_index) ./ expA_set(start_index:end_index); % norm
simA_speed_step = simA_speed(start_index:end_index) ./ simA_set(start_index:end_index); % norm

% make ideal response data for tau = 36 ms
theoretical_time = linspace(0, time_limit, 1000);
theoretical_response = 1 - exp(-theoretical_time / tau); % ideal step response

% key times for visualisation
key_times = [tau, 2*tau, 3*tau, 4*tau, 5*tau];
key_values = 1 - exp(-key_times / tau); % ideal values at key times

% Plot
figure;
plot(timeA_step, expA_speed_step, 'b', 'LineWidth', 1.5); hold on; % Exp response
plot(timeA_step, simA_speed_step, 'g--', 'LineWidth', 1.5); % Sim response
plot(theoretical_time, theoretical_response_offset, 'r--', 'LineWidth', 1.5); % Theory response

% key points on plot
plot(key_times, key_values, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5); % Key points
for i = 1:length(key_times)
    text(key_times(i), key_values(i), ...
        ['\leftarrow ', num2str(key_times(i)), ' s'], 'FontSize', 10, 'Color', 'k');
end

% Add labels, title, and legend
xlabel('Time Since Step (s)');
ylabel('Normalised Speed/Set-Speed Response');
title('Response Plot for (0.3 s After Step)');
legend('Experimental Response', 'Simulated Response', ...
       'Theoretical Response', 'Key Points (\tau, 2\tau, 3\tau, 3\tau, 4\tau, 5\tau)');
grid on;