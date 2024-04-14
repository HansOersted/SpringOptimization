close all

% step 1. generate 50 random points within a unit circle

% Number of points
N = 50;

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

% Plot the points
figure;
plot(x, y, 'o');
axis equal;
grid on;
xlabel('X coordinate');
ylabel('Y coordinate');
title('50 Random Points Inside a Unit Circle Centered at (0, -1)');

% Add a circle to verify the boundary
hold on;
theta = linspace(0, 2*pi, 100);
circle_x = cos(theta);
circle_y = sin(theta) - 1; % Adjusting y coordinates for the unit circle
plot(circle_x, circle_y, 'r--'); % Plotting unit circle boundary
hold off;

% step 2. Generate triangles of delaunay and determine the springs 

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

% step 3. 准备好仿真器
% step 4. 计算平衡后的情况（可以利用仿真器）
% step 5. 评估平衡后的情况，比如方不方
% step 6. 调整所有的弹簧的k与L0，重复step3 -- 5