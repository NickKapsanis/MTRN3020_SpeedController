
function [A, tau_m] = motorData() 
    % first find the motor gain A, and the motor time constant t_m
    run('C:\Users\nickk\OneDrive\Desktop\MTRN3020\LAB2_SpeedController\Laboratory Instructions-20250414\noload.m');

    time = TEST(:,1) / 1000; % time vector (ms) -> seconds
    speed = TEST(:,3); % speed vector (Counts/s)

    %plot the speed response
    figure;
    plot(time, speed, 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Speed (Counts/s)');
    title('Motor Speed Response');
    grid on;

    % we have the transfer function of the motor in the s-domain, 
    % we want the step response of the motor in the s-domain 
    model = @(params, t) params(1) * (1 - exp(-t / params(2))); % params(1) = A, params(2) = tau_m
    initial_guess = [max(speed), 0.1]; % Guess A as max speed, tau_m as 0.1s
    params = lsqcurvefit(model, initial_guess, time, speed);
    A = params(1);
    tau_m = params(2);

    A = A / 24; %A is calculated here with a unit step response, but we actually have a 24V step response, so we need to divide by 24V

    % Display results
    fprintf('Motor Gain (A): %.4f counts/s per volt\n', A);
    fprintf('Motor Time Constant (tau_m): %.4f s\n', tau_m);

    % Plot fitted curve
    hold on;
    fitted_speed = model(params, time);
    plot(time, fitted_speed, 'r--', 'LineWidth', 1.5);
    legend('Measured Speed', 'Fitted Model');
end