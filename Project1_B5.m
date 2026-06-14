% B.5: Compute kinematic quantities over full cycle

% Define known variables and equations from Part B
%   as θ3, θ4, etc... will need to be solved for
%% Links
r1 = 8.0; %cm
r2 = 3.0; %cm
r3 = 10.0; %cm
r4 = 9.0; %cm

w2 = 36; %rad/sec
a2 = 0; %rad/sec^2

rbp = 5.0; %cm
r5 = rbp;

beta = 30; %degrees
beta_rad = deg2rad(beta);

%% H values
h1 = r1/r2;
h2 = r1/r3;
h3 = r1/r4;
h4 = (-r1^2-r2^2-r3^2+r4^2)/(2*r2*r3);
h5 = (r1^2+r2^2-r3^2+r4^2)/(2*r2*r4);

theta2 = 0:2:360; %input incremements of 2 degrees over a full rotation
theta2_rad = deg2rad(theta2);
theta2dot_rad = 36;
theta2doubledot_rad = 0;

%% Constants for simplification
a = h4 + h2.*cos(theta2_rad)-h1+cos(theta2_rad);
b = -sin(theta2_rad);
c = h4+h2.*cos(theta2_rad)+h1-cos(theta2_rad);
d = h5-h3.*cos(theta2_rad)+cos(theta2_rad)-h1;
e = h5-h3.*cos(theta2_rad)-cos(theta2_rad)+h1;

%% Desired Angle for Verification
indexAt60 = find(theta2 == 60);

%% a) θ3, θ4 vs θ2
theta3 = 2.*atan((-b-sqrt(b.^2-a.*c))./(a));
theta3_deg = rad2deg(theta3);

theta4 = 2.*atan((-b-sqrt(b.^2-d.*e))./(d));
theta4_deg = rad2deg(theta4);

theta3_deg_60 = theta3_deg(indexAt60);
theta4_deg_60 = theta4_deg(indexAt60);

figure ('Name', 'Ang.Disp.');
plot(theta2,theta3_deg);
title('Angular Displacements in degrees');
xlabel('θ2 [°]');
ylabel('θ3 & θ4 [°]');
hold on;
plot(theta2,theta4_deg);

plot(60, theta3_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,theta3_deg_60 + 1, sprintf('  (%.1f, %.3f)', 60,theta3_deg_60));
plot(60, theta4_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,theta4_deg_60, sprintf('  (%.1f, %.3f)', 60,theta4_deg_60));

legend({'θ3','θ4'},'Location','Southeast');
hold off;

%% b) θ̇3, θ̇4 vs θ2
theta3_dot = (r2.*theta2dot_rad.*sin(theta2_rad-theta4))./ ...
    (r3.*sin(theta4-theta3));
theta4_dot = (r2.*theta2dot_rad.*sin(theta2_rad-theta3))./ ...
    (r4.*sin(theta4-theta3));

theta3_dot_deg_60 = theta3_dot(indexAt60);
theta4_dot_deg_60 = theta4_dot(indexAt60);

figure ('Name', 'Ang.Vel.');
plot(theta2,theta3_dot);
title('Angular Velocities in rad/sec');
xlabel('θ2 [°]');
ylabel('θ̇3 & θ̇4 [rad/sec]');
hold on;
plot(theta2,theta4_dot);

plot(60, theta3_dot_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,theta3_dot_deg_60, sprintf('  (%.1f, %.3f)', 60,theta3_dot_deg_60));
plot(60, theta4_dot_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,theta4_dot_deg_60, sprintf('  (%.1f, %.3f)', 60,theta4_dot_deg_60));

legend({'θ̇3','θ̇4'},'Location','Southeast');
hold off;

%% c) θ̈3, θ̈4 vs θ2
theta3_doubledot = ( ...
    -r4.*theta4_dot.^2 ...
    + r3.*theta3_dot.^2.*cos(theta3 - theta4) ...
    + r2.*theta2dot_rad.^2.*cos(theta2_rad - theta4) ...
    + r2.*theta2doubledot_rad.*sin(theta2_rad - theta4) ...
    ) ./ (r3.*sin(theta4 - theta3));
theta4_doubledot = (-r3.*theta3_dot.^2 ...
    + r4.*theta4_dot.^2.*cos(theta4 - theta3) ...
    - r2.*theta2dot_rad.^2.*cos(theta2_rad - theta3) ...
    + r2.*theta2doubledot_rad.*sin(theta2_rad - theta3) ...
    ) ./ (r4.*sin(theta3 - theta4));

theta3_doubledot_deg_60 = theta3_doubledot(indexAt60);
theta4_doubledot_deg_60 = theta4_doubledot(indexAt60);

figure ('Name', 'Ang.Acc.');
plot(theta2,theta3_doubledot);
title('Angular Accelerations in rad/sec²');
xlabel('θ2 [°]');
ylabel('θ̈3 & θ̈4 [rad/sec²]');
hold on;
plot(theta2,theta4_doubledot);

plot(60, theta3_doubledot_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,theta3_doubledot_deg_60, sprintf('  (%.1f, %.3f)', 60,theta3_doubledot_deg_60));
plot(60, theta4_doubledot_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,theta4_doubledot_deg_60, sprintf('  (%.1f, %.3f)', 60,theta4_doubledot_deg_60));

legend({'θ̈3','θ̈4'},'Location','Southeast');
hold off;

%% d) xp, yp vs θ2
theta5 = theta3+beta_rad;

xp = r5.*cos(theta5) + r2.*cos(theta2_rad);
yp = r5.*sin(theta5) + r2.*sin(theta2_rad);

xp_deg_60 = xp(indexAt60);
yp_deg_60 = yp(indexAt60);

figure ('Name', 'Cart.Pos.');
plot(theta2,xp);
title('Cartesian Position Components in cm');
xlabel('θ2 [°]');
ylabel('xp & yp [cm]');
hold on;
plot(theta2,yp);

plot(60, xp_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,xp_deg_60, sprintf('  (%.1f, %.3f)', 60,xp_deg_60));
plot(60, yp_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,yp_deg_60, sprintf('  (%.1f, %.3f)', 60,yp_deg_60));

legend({'xp','yp'},'Location','Southeast');
hold off;

%% e) Velocity magnitude |vp| vs θ2
xp_dot = -theta3_dot.*r5.*sin(theta5) ...
         -theta2dot_rad.*r2.*sin(theta2_rad);
yp_dot = theta3_dot.*r5.*cos(theta5) ...
         +theta2dot_rad.*r2.*cos(theta2_rad);

xp_dot_deg_60 = xp_dot(indexAt60);
yp_dot_deg_60 = yp_dot(indexAt60);
vp_mag_deg_60 = sqrt(xp_dot_deg_60.^2+yp_dot_deg_60.^2);

vp_mag = abs(sqrt(xp_dot.^2 + yp_dot.^2));
figure ('Name', 'Vel.Mag.');
plot(theta2,vp_mag);
title('Velocity Magnitude in cm/sec²');
xlabel('θ2 [°]');
ylabel('vp [cm]');
hold on;

plot(60, vp_mag_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,vp_mag_deg_60, sprintf('  (%.1f, %.3f)', 60,vp_mag_deg_60));

legend({'vp'},'Location','Southeast');
hold off;

%% f) Acceleration Magnitude |ap| vs θ2
xp_doubledot = -theta3_doubledot.*r5.*sin(theta5) ...
               -theta3_dot.^2.*r5.*cos(theta5) ...
               -theta2doubledot_rad.*r2.*sin(theta2_rad) ...
               -theta2dot_rad.^2.*r2.*cos(theta2_rad);
yp_doubledot = theta3_doubledot.*r5.*cos(theta5) ...
               -theta3_dot.^2.*r5.*sin(theta5) ...
               +theta2doubledot_rad.*r2.*cos(theta2_rad) ...
               -theta2dot_rad.^2.*r2.*sin(theta2_rad);

xp_doubledot_deg_60 = xp_doubledot(indexAt60);
yp_doubledot_deg_60 = yp_doubledot(indexAt60);
ap_mag_deg_60 = sqrt(xp_doubledot_deg_60.^2+yp_doubledot_deg_60.^2);

ap_mag = abs(sqrt(xp_doubledot.^2 + yp_doubledot.^2));
figure ('Name', 'Acc.Mag.');
plot(theta2, ap_mag);
title('Acceleration Magnitude in cm/sec²');
xlabel('θ2 [°]');
ylabel('ap [cm]');
hold on;

plot(60, ap_mag_deg_60, 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
text(60,ap_mag_deg_60, sprintf('  (%.1f, %.3f)', 60,ap_mag_deg_60));

legend({'ap'}, 'Location', 'Southeast');
hold off;

%% g) Coupler Curve yp vs xp
figure ('Name', 'Coupler Curve');
plot(xp, yp);
title('Coupler Curve: yp vs xp');
xlabel('xp [cm]');
ylabel('yp [cm]');
grid on;