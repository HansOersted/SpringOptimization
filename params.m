close all
clear

%% Define the parameters of the spring and the particle

g = 9.81;
m = 1.0;  % mass of the particle

L0 = 1.0;  % initial length of the spring (nature length) 
k = 10;  % spring constant
beta = 2;  % damping constant


initial_theta = 0.5;
initial_L = 0;

simOut = sim(spring_simulator, 'SimulationMode', 'normal', 'SrcWorkspace', 'current');

%% Play the video

% Prepare the data from simulink
time = simOut.position.time;
position_values = simOut.position.signals.values;

% receive x and y from position_values
x = position_values(:, 1);
y = position_values(:, 2);

ref_x = 0.0;
ref_y = 0.0;

h_fig = figure;
hold on;
h_plot = plot(x(1), y(1), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
h_ref = plot(ref_x, ref_y, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k'); % reference point

% line linking the points
h_line1 = line([x(1), ref_x], [y(1), ref_y], 'Color', 'r', 'LineStyle', '-');

xlabel('X Position');
ylabel('Y Position');
title('Spring Mass System');
legend([h_ref, h_plot], 'm0', 'm1');
axis equal;
grid on;

% axix
x_all = [x; ref_x];
y_all = [y; ref_y];
xlim([min(x_all), max(x_all)]);
ylim([min(y_all), max(y_all)]);

% play the video
for k = 1:length(time)
    % updata the position
    set(h_plot, 'XData', x(k), 'YData', y(k));
    
    % update the line
    set(h_line1, 'XData', [x(k), ref_x], 'YData', [y(k), ref_y]);

    % fresh and pause
    drawnow;
    pause(0.00001); 
end