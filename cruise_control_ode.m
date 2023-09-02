function y = cruise_control_ode(t, V, desired_speed, gain, time_constant)
%     time_constant = 100; % Time constant of the vehicle in seconds
%     gain = 20; % Gain of the controller

    % Define the PI controller gains
    Kp = 50;
%     Kp = 0.4;  % Proportional gain
    Ki = 0.02;  % Integral gain
    
    % Define the PI controller function
    error_integral = 0;
    Throttle = 0.3 + Kp*(desired_speed - V) + Ki*error_integral;
%     pi_controller = @(V) 0.3 + Kp*(desired_speed - V) + Ki*error_integral;
%     
%     Throttle = pi_controller(V);

    drag = 0.01; % Drag coefficient of the vehicle

    dVdt = (1/time_constant) * (Throttle - gain*V - drag*V'*V) + (gain/time_constant)*(desired_speed - V);
    y = dVdt;
end
