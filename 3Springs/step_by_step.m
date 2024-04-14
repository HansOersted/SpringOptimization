close all
clear
%% step 1. generate 3 random points within a unit circle

% Number of points
N = 3;

% Path to save or load the data
dataPath = 'random_points.mat';

% Check if the file exists
if exist(dataPath, 'file')
    % Load saved points
    load(dataPath, 'x', 'y');
else
    % Initialize arrays for Cartesian coordinates
    x = zeros(N, 1);
    y = zeros(N, 1);

    % Generate points
    for i = 1:N
        r = sqrt(rand());  % Generate radius (square root ensures uniform distribution)
        theta = 2 * pi * rand();  % Generate angle

        % Convert polar coordinates to Cartesian coordinates and adjust y coordinate
        x(i) = r * cos(theta);
        y(i) = r * sin(theta) - 1; % Adjusting y coordinates to shift the circle center to (0, -1)
    end
    
    % Save the generated points
    save(dataPath, 'x', 'y');
end

% Plot the points
figure;
plot(x, y, 'o');
axis equal;
grid on;
xlabel('X coordinate');
ylabel('Y coordinate');
title('3 Random Points Inside a Unit Circle Centered at (0, -1)');

% Add a circle to verify the boundary
hold on;
theta = linspace(0, 2*pi, 100);
circle_x = cos(theta);
circle_y = sin(theta) - 1; % Adjusting y coordinates for the unit circle
plot(circle_x, circle_y, 'r--'); % Plotting unit circle boundary
hold off;

%% step 2. Perform Delaunay triangulation

% Perform Delaunay triangulation
tri = delaunay(x, y);

% Plot the triangulation
figure;
triplot(tri, x, y);
axis equal;
grid on;
xlabel('X coordinate');
ylabel('Y coordinate');
title('Delaunay Triangulation of Points Inside a Unit Circle at (0, -1)');


%% step 3. Initial conditions for the simulator

% Sort the points by y-coordinate in descending order
[sorted_y, idx] = sort(y, 'descend'); % 'descend' ensures sorting is from highest to lowest
sorted_x = x(idx);

% hold the coordinates of the point with the highest y-coordinate,
ref_x = sorted_x(1); 
ref_y = sorted_y(1);

% initial_x1, initial_y1 will hold the coordinates of the second highest, and so on.
initial_x1 = sorted_x(2);
initial_y1 = sorted_y(2);
initial_x2 = sorted_x(3);
initial_y2 = sorted_y(3);

Lr1 = sqrt((initial_x1 - ref_x)^2 + (initial_y1 - ref_y)^2);
Lr2 = sqrt((initial_x2 - ref_x)^2 + (initial_y2 - ref_y)^2);
Lr3 = sqrt((initial_x1 - initial_x2)^2 + (initial_y1 - initial_y2)^2);

m1 = 1;
m2 = 1;

g = 10;

% waiting to be updated ...
k1 = 10;
k2 = 40;
k3 = 10;

% waiting to be updated ...
beta1 = 1;
beta2 = 1;
beta3 = 1;

simOut = sim(simulator, 'SimulationMode', 'normal', 'SrcWorkspace', 'current');

%% step 4. Play the video

% Prepare the data from simulink
time = simOut.positions.time;
position_values = simOut.positions.signals.values;

% receive x and y from position_values
x1 = position_values(:, 1);
y1 = position_values(:, 2);
x2 = position_values(:, 3);
y2 = position_values(:, 4);

h_fig = figure;
hold on;
h_plot1 = plot(x1(1), y1(1), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
h_plot2 = plot(x2(1), y2(1), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
h_ref = plot(ref_x, ref_y, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k'); % reference point

% line linking the points
h_line1 = line([x1(1), ref_x], [y1(1), ref_y], 'Color', 'b', 'LineStyle', '-');
h_line2 = line([x2(1), ref_x], [y2(1), ref_y], 'Color', 'r', 'LineStyle', '-');
h_line3 = line([x1(1), x2(1)], [y1(1), y2(1)], 'Color', 'k', 'LineStyle', '-');

xlabel('X Position');
ylabel('Y Position');
title('2D Position Animation with Connecting Lines');
legend([h_ref, h_plot1, h_plot2], 'm0', 'm1', 'm2');
axis equal;
grid on;

% axis
x_all = [x1; x2; ref_x];
y_all = [y1; y2; ref_y];
xlim([min(x_all), max(x_all)]);
ylim([min(y_all), max(y_all)]);

% play the video
for k = 1:length(time)
    % updata the position
    set(h_plot1, 'XData', x1(k), 'YData', y1(k));
    set(h_plot2, 'XData', x2(k), 'YData', y2(k));
    
    % update the line
    set(h_line1, 'XData', [x1(k), ref_x], 'YData', [y1(k), ref_y]);
    set(h_line2, 'XData', [x2(k), ref_x], 'YData', [y2(k), ref_y]);
    set(h_line3, 'XData', [x1(k), x2(k)], 'YData', [y1(k), y2(k)]);

    % fresh and pause
    drawnow;
    pause(0.00001); 
end


%% step 5. Evaluate the triangle

x1_final = x1(end);
y1_final = y1(end);

x2_final = x2(end);
y2_final = y2(end);

LL1 = sqrt( (x1_final-ref_x)^2 + (y1_final-ref_y)^2 );

LL2 = sqrt( (x2_final-ref_x)^2 + (y2_final-ref_y)^2 );

LL3 = sqrt( (x1_final-x2_final)^2 + (y1_final-y2_final)^2 );

mean = (LL1 + LL2 + LL3) / 3;

cost = (LL1 - mean)^2 + (LL2 - mean)^2 + (LL3 - mean)^2;

%% step 6. Change k1 , k2, k3 and repeat step 3 -- 5



