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
grid on;
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
grid on;
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
grid on;
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
grid on;
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
grid on;
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
grid on;
hold off;

%% g) Coupler Curve yp vs xp
figure ('Name', 'Coupler Curve');
plot(xp, yp);
title('Coupler Curve: yp vs xp');
xlabel('xp [cm]');
ylabel('yp [cm]');
grid on;

%% C.3: Complete Force Analysis over one crank revolution
m2 = 8; %g
m3 = 30; %g
m4 = 20; %g

b2 = 1.50; %cm to mass center
b3 = 4.00; %cm to mass center
b4 = 3.50; %cm to mass center

I_G3 = 250; %gr*cm^2
I_G4 = 90; %gr*cm^2

%% Positions and Accelerations of mass centers
% Link 2: Position
X_G2 = b2.*cos(theta2_rad);
Y_G2 = b2.*sin(theta2_rad);

% Link 2: Acceleration
A_X_G2 = -b2.*theta2doubledot_rad.*sin(theta2_rad)-b2.*theta2dot_rad.^2 ...
    .*cos(theta2_rad);
A_Y_G2 = b2.*theta2doubledot_rad.*cos(theta2_rad)-b2.*theta2dot_rad.^2 ...
    .*sin(theta2_rad);

% Link 3: Position
X_G3 = r2.*cos(theta2_rad) + b3.*cos(theta3);
Y_G3 = r2.*sin(theta2_rad) + b3.*sin(theta3);

% Link 3: Accceleration
A_X_G3 = -r2.*theta2doubledot_rad.*sin(theta2_rad) ...
    -r2.*theta2dot_rad.^2.*cos(theta2_rad) ...
    -b3.*theta3_doubledot.*sin(theta3) ...
    -b3.*theta3_dot.^2.*cos(theta3);
A_Y_G3 = r2.*theta2doubledot_rad.*cos(theta2_rad) ...
    -r2.*theta2dot_rad.^2.*sin(theta2_rad) ...
    +b3.*theta3_doubledot.*cos(theta3) ...
    -b3.*theta3_dot.^2.*sin(theta3);

% Link 4: Position
X_G4 = r1 + b4.*cos(theta4);
Y_G4 = b4.*sin(theta4);

% Link 4: Acceleration
A_X_G4 = -b4.*theta4_doubledot.*sin(theta4)-b4.*theta4_dot.^2 ...
    .*cos(theta4);
A_Y_G4 = b4.*theta4_doubledot.*cos(theta4)-b4.*theta4_dot.^2 ...
    .*sin(theta4);

n = length(theta2_rad);
N_to_cm_g = 1e5; %g*cm/s^2
Fext = 0.5*N_to_cm_g; %N

F12x = zeros(1,n);
F12y = zeros(1,n);
F32x = zeros(1,n);
F32y = zeros(1,n);
F34x = zeros(1,n);
F34y = zeros(1,n);
F41x = zeros(1,n);
F41y = zeros(1,n);
M12 = zeros(1,n); % Moment applied at 02, same as W in notes

I_O4 = I_G4 + m4*b4^2; % parallel axis theorem moved from I_G4

%% Calculate and fill forces & moment from matrix
for k = 1:n
    A = [ ...
        -1,  0,  1,  0,  0,  0,  0,  0,  0;
         0, -1,  0,  1,  0,  0,  0,  0,  0;
         r2*sin(theta2_rad(k)), -r2*cos(theta2_rad(k)), 0, 0, 0, 0, 0, 0, 1;
         0,  0, -1,  0, -1,  0,  0,  0,  0;
         0,  0,  0, -1,  0, -1,  0,  0,  0;
         0,  0, -b3*sin(theta3(k)), b3*cos(theta3(k)), ...
                 b3*sin(theta3(k)), -b3*cos(theta3(k)), 0, 0, 0;
         0,  0,  0,  0,  1,  0, -1,  0,  0;
         0,  0,  0,  0,  0,  1,  0, -1,  0;
         0,  0,  0,  0, -r4*sin(theta4(k)), r4*cos(theta4(k)), 0, 0, 0 ...
        ];
    B = [ ...
        m2*A_X_G2(k);
        m2*A_Y_G2(k);
        0;
        m3*A_X_G3(k);
        m3*A_Y_G3(k) + Fext;
        I_G3*theta3_doubledot(k) - Fext*(b3*cos(theta3(k)) - rbp*cos(theta3(k)));
        m4*A_X_G4(k);
        m4*A_Y_G4(k);
        I_O4*theta4_doubledot(k) ...
        ];
    x = A\B; % solve for temporary variable

    F12x(k) = x(1);
    F12y(k) = x(2);
    F32x(k) = x(3);
    F32y(k) = x(4);
    F34x(k) = x(5);
    F34y(k) = x(6);
    F41x(k) = x(7);
    F41y(k) = x(8);
    M12(k) = x(9);
end

%% (i) Driving Torque M12 vs θ2
figure ('Name', 'Driving Torque');
M12_N = M12./N_to_cm_g;
plot(theta2, M12_N);
title('Driving Torque: M12 vs θ2');
xlabel('θ2 [°]');
ylabel('M12 [Nm]');
grid on;

%% (ii) Polar Plot of shaking force Fs (mag and dir)
Fs_mag = sqrt((-F12x+F41x).^2+(-F12y+F41y).^2)./N_to_cm_g; %txtbook p.391 9.2-12
Fs_dir = atan2((-F12y+F41y), (-F12x+F41x)); %txtbook p.391 9.2-10

figure('Name', 'Polar Fs');
polarplot(Fs_dir, Fs_mag);
title('Shaking Force [N] at θ2');
legend({'Fs [N]'})
grid on;

%% (iii) Components of Fs vs θ2
Fs_x = (-F12x+F41x)./N_to_cm_g;
Fs_y = (-F12y+F41y)./N_to_cm_g;

figure('Name', 'Shaking Force Componenets');
plot(theta2, Fs_x, 'b', theta2, Fs_y, 'r');
title('Shaking Force Components vs θ2');
xlabel('θ2 [°]');
ylabel('Force Components [N]');
legend({'Fs_x [N]', 'Fs_y [N]'}, 'Location','southeast');
grid on;
hold off;

%% (iv) Shaking Moment Ms vs θ2
Ms = M12_N + ((F41y).*r1./N_to_cm_g);

figure('Name','Shaking Moment');
plot(theta2,Ms);
title('Shaking Moment vs θ2');
xlabel('θ2 [°]');
ylabel('Ms [Nm]');
grid on;

%% (v) Force Magnitude at each turning pair vs θ2
F12_mag = sqrt(F12x.^2 + F12y.^2)./N_to_cm_g;
F23_mag = sqrt(F32x.^2 + F32y.^2)./N_to_cm_g;
F34_mag = sqrt(F34x.^2 + F34y.^2)./N_to_cm_g;
F14_mag = sqrt(F41x.^2 + F41y.^2)./N_to_cm_g;

figure('Name','Turning Pairs');
plot(theta2,F12_mag,'b',theta2,F23_mag,'g',theta2,F34_mag,'r',theta2,F14_mag,'y');
title('Turning Pair Force Magnitude vs θ2');
xlabel('θ2 [°]');
ylabel('Turning Pairs [N]');
legend({'|F12| [N]','|F23| [N]','|F34| [N]','|F14| [N]'},'Location','northeast');
grid on;
