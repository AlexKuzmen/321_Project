%% MTE 321 Project 1 - Four-Bar Mechanism Analysis
% Version 5.0: added Part D) Modularized Code, added keywords to replace
% symbols (ex: θ2 > \theta_2, 'latex-like')
% Units used:
%   - length: cm
%   - mass: g
%   - time: s
%   - force internal matrix units: g*cm/s^2
%   - plotted force outputs: N
%   - plotted moment outputs: N*cm
 clear vars; close all; clc;
%% SETUP VARIABLES
% Define all constant parameters to given and to be used
p = defineParameters(); % changed to p from parameters to improve read

%% B.5 (a), (b), (c) using Loop 1 Equations (see report)
pos = angularDisplacements(p); % θ3, θ4
vel = angularVelocities(p, pos); % θ̇3, θ̇4
acc = angularAccelerations(p, pos, vel); % θ̈3, θ̈4

%% B.5 (d), (e), (f), (g) using Loop 2 Equations (see report)
pointP = couplerPointAnalysis(p, pos, vel, acc);

% Generate Plots for B
plotPartB(p, pos, vel, acc, pointP); 

%% C.3 force analysis with external force applied (see report)
forceLoaded = forceAnalysis(p, pos, vel, acc, p.Fext_N);

% Generate Plots for C
plotPartC(p, forceLoaded);

%% C.4 no loaded force applied
forceNoLoad = forceAnalysis(p, pos, vel, acc, 0);

% Print numerical C.4 answers to the MATLAB Command Window.
C4 = summarizeC4(p, forceLoaded, forceNoLoad);
printC4Summary(C4);

%%__________________________Functions______________________________________
function p = defineParameters()
% Link lengths [cm]
    p.r1 = 8.0;
    p.r2 = 3.0;
    p.r3 = 10.0;
    p.r4 = 9.0;

    % Input crank motion
    p.w2 = 36; %rad/s
    p.a2 = 0; %rad/s^2

    % Coupler point P location relative to point B
    p.rbp = 5.0; %cm
    p.beta_deg = 30; %deg
    p.beta_rad = deg2rad(p.beta_deg); %rad

    % Masses
    p.m2 = 8; %g
    p.m3 = 30; %g
    p.m4 = 20; %g

    p.b2 = 1.50; %cm from O2 to G2
    p.b3 = 4.00; %cm from B to G3
    p.b4 = 3.50; %cm from O4 to G4

    p.I_G3 = 250; %g*cm^2
    p.I_G4 = 90; %g*cm^2

    % External load at point P
    p.Fext_N = 0.5; %N

    % Unit conversions
    p.N_to_gcm = 1e5; %1 N = 1e5 g*cm/s^2
    p.Ncm_to_Nm = 1/100; %1 N*cm = 0.01 N*m

    % Input crank angle range for one full cycle
    p.theta2_deg = 0:2:360; %increments of 2 deg
    p.theta2_rad = deg2rad(p.theta2_deg); %rad
    p.indexAt60 = find(p.theta2_deg == 60, 1); %sanity check for 1 match
end

function pos = angularDisplacements(p)
    %% B.5(a): Loop 1 angular displacements
    % H Constants
    h1 = p.r1/p.r2;
    h2 = p.r1/p.r3;
    h3 = p.r1/p.r4;
    h4 = (-p.r1^2 - p.r2^2 - p.r3^2 + p.r4^2)/(2*p.r2*p.r3);
    h5 = ( p.r1^2 + p.r2^2 - p.r3^2 + p.r4^2)/(2*p.r2*p.r4);

    % Quadratic substitution simplifiers
    a = h4 + h2.*cos(p.theta2_rad) - h1 + cos(p.theta2_rad);
    b = -sin(p.theta2_rad);
    c = h4 + h2.*cos(p.theta2_rad) + h1 - cos(p.theta2_rad);
    d = h5 - h3.*cos(p.theta2_rad) + cos(p.theta2_rad) - h1;
    e = h5 - h3.*cos(p.theta2_rad) - cos(p.theta2_rad) + h1;

    % Simplified & cancelled quadratic (see report)
    pos.theta3 = 2.*atan((-b - sqrt(b.^2 - a.*c))./a);
    pos.theta4 = 2.*atan((-b - sqrt(b.^2 - d.*e))./d);

    pos.theta3_deg = rad2deg(pos.theta3);
    pos.theta4_deg = rad2deg(pos.theta4);
end

function vel = angularVelocities(p, pos)
    %% B.5(b): Loop 1 angular velocities
    vel.theta3_dot = (p.r2.*p.w2.*sin(p.theta2_rad - pos.theta4))./ ...
        (p.r3.*sin(pos.theta4 - pos.theta3));

    vel.theta4_dot = (p.r2.*p.w2.*sin(p.theta2_rad - pos.theta3))./ ...
        (p.r4.*sin(pos.theta4 - pos.theta3));
end

function acc = angularAccelerations(p, pos, vel)
    %% B.5(c): Loop 1 angular accelerations
    acc.theta3_doubledot = ( ...
        -p.r4.*vel.theta4_dot.^2 ...
        + p.r3.*vel.theta3_dot.^2.*cos(pos.theta3 - pos.theta4) ...
        + p.r2.*p.w2.^2.*cos(p.theta2_rad - pos.theta4) ...
        + p.r2.*p.a2.*sin(p.theta2_rad - pos.theta4) ...
        ) ./ (p.r3.*sin(pos.theta4 - pos.theta3));

    acc.theta4_doubledot = ( ...
        -p.r3.*vel.theta3_dot.^2 ...
        + p.r4.*vel.theta4_dot.^2.*cos(pos.theta4 - pos.theta3) ...
        - p.r2.*p.w2.^2.*cos(p.theta2_rad - pos.theta3) ...
        + p.r2.*p.a2.*sin(p.theta2_rad - pos.theta3) ...
        ) ./ (p.r4.*sin(pos.theta3 - pos.theta4));
end

function pointP = couplerPointAnalysis(p, pos, vel, acc)
    %% B.5(d) to (g): Coupler point P position, velocity, acceleration
    pointP.theta5 = pos.theta3 + p.beta_rad;

    % Position Components of P [cm]
    pointP.xp = p.rbp.*cos(pointP.theta5) + p.r2.*cos(p.theta2_rad);
    pointP.yp = p.rbp.*sin(pointP.theta5) + p.r2.*sin(p.theta2_rad);

    % Velocity components of P [cm/s]
    pointP.xp_dot = -vel.theta3_dot.*p.rbp.*sin(pointP.theta5) ...
                -p.w2.*p.r2.*sin(p.theta2_rad);
    pointP.yp_dot =  vel.theta3_dot.*p.rbp.*cos(pointP.theta5) ...
                +p.w2.*p.r2.*cos(p.theta2_rad);
    pointP.vp_mag = sqrt(pointP.xp_dot.^2 + pointP.yp_dot.^2);

    % Acceleration components of P [cm/s^2]
    pointP.xp_doubledot = ...
        -acc.theta3_doubledot.*p.rbp.*sin(pointP.theta5) ...
                      -vel.theta3_dot.^2.*p.rbp.*cos(pointP.theta5) ...
                      -p.a2.*p.r2.*sin(p.theta2_rad) ...
                      -p.w2.^2.*p.r2.*cos(p.theta2_rad);

    pointP.yp_doubledot =  ...
        acc.theta3_doubledot.*p.rbp.*cos(pointP.theta5) ...
                      -vel.theta3_dot.^2.*p.rbp.*sin(pointP.theta5) ...
                      +p.a2.*p.r2.*cos(p.theta2_rad) ...
                      -p.w2.^2.*p.r2.*sin(p.theta2_rad);

    % Magnitude of acceleration of point P
    pointP.ap_mag = sqrt(pointP.xp_doubledot.^2 + pointP.yp_doubledot.^2);
end

function force = forceAnalysis(p, pos, vel, acc, Fext_N)
    %% Part C.3: Complete force analysis over one full crank revolution
    n = length(p.theta2_rad); % define vector forces to fill
    Fext = Fext_N*p.N_to_gcm; % convert N to g*cm/s^2

    % Link 2 mass-center position and acceleration
    force.X_G2 = p.b2.*cos(p.theta2_rad);
    force.Y_G2 = p.b2.*sin(p.theta2_rad);
    force.A_X_G2 = -p.b2.*p.a2.*sin(p.theta2_rad) ...
                   -p.b2.*p.w2.^2.*cos(p.theta2_rad);
    force.A_Y_G2 =  p.b2.*p.a2.*cos(p.theta2_rad) ...
                   -p.b2.*p.w2.^2.*sin(p.theta2_rad);

    % Link 3 mass-center position and acceleration
    force.X_G3 = p.r2.*cos(p.theta2_rad) + p.b3.*cos(pos.theta3);
    force.Y_G3 = p.r2.*sin(p.theta2_rad) + p.b3.*sin(pos.theta3);
    force.A_X_G3 = -p.r2.*p.a2.*sin(p.theta2_rad) ...
                   -p.r2.*p.w2.^2.*cos(p.theta2_rad) ...
                   -p.b3.*acc.theta3_doubledot.*sin(pos.theta3) ...
                   -p.b3.*vel.theta3_dot.^2.*cos(pos.theta3);
    force.A_Y_G3 =  p.r2.*p.a2.*cos(p.theta2_rad) ...
                   -p.r2.*p.w2.^2.*sin(p.theta2_rad) ...
                   +p.b3.*acc.theta3_doubledot.*cos(pos.theta3) ...
                   -p.b3.*vel.theta3_dot.^2.*sin(pos.theta3);

    % Link 4 mass-center position and acceleration
    force.X_G4 = p.r1 + p.b4.*cos(pos.theta4);
    force.Y_G4 = p.b4.*sin(pos.theta4);
    force.A_X_G4 = -p.b4.*acc.theta4_doubledot.*sin(pos.theta4) ...
                   -p.b4.*vel.theta4_dot.^2.*cos(pos.theta4);
    force.A_Y_G4 =  p.b4.*acc.theta4_doubledot.*cos(pos.theta4) ...
                   -p.b4.*vel.theta4_dot.^2.*sin(pos.theta4);

    % Initialize force and applied torque with zeroes
    force.F12x = zeros(1,n);
    force.F12y = zeros(1,n);
    force.F32x = zeros(1,n);
    force.F32y = zeros(1,n);
    force.F34x = zeros(1,n);
    force.F34y = zeros(1,n);
    force.F41x = zeros(1,n);
    force.F41y = zeros(1,n);
    force.M12  = zeros(1,n);

    % Parallel-axis theorem for link 4 inertia about O4
    I_O4 = p.I_G4 + p.m4*p.b4^2;

    % Solve the 9 x 9 force system at every theta2 value
    % For a cleaner matrix, see report
    for k = 1:n
        A = [ ...
            -1,  0,  1,  0,  0,  0,  0,  0,  0;
             0, -1,  0,  1,  0,  0,  0,  0,  0;
             p.r2*sin(p.theta2_rad(k)), -p.r2*cos(p.theta2_rad(k)), ...
             0, 0, 0, 0, 0, 0, 1;
             0,  0, -1,  0, -1,  0,  0,  0,  0;
             0,  0,  0, -1,  0, -1,  0,  0,  0;
             0,  0, -p.b3*sin(pos.theta3(k)), p.b3*cos(pos.theta3(k)), ...
                     p.b3*sin(pos.theta3(k)), -p.b3*cos(pos.theta3(k)), ...
                     0, 0, 0;
             0,  0,  0,  0,  1,  0, -1,  0,  0;
             0,  0,  0,  0,  0,  1,  0, -1,  0;
             0,  0,  0,  0, -p.r4*sin(pos.theta4(k)), ...
             p.r4*cos(pos.theta4(k)), 0, 0, 0 ...
            ];

        % Known Vecotr B (see report) 
        B = [ ...
            p.m2*force.A_X_G2(k);
            p.m2*force.A_Y_G2(k);
            0;
            p.m3*force.A_X_G3(k);
            p.m3*force.A_Y_G3(k) + Fext;
            p.I_G3*acc.theta3_doubledot(k) - ...
            Fext*(p.b3*cos(pos.theta3(k)) - p.rbp*cos(pos.theta3(k)));
            p.m4*force.A_X_G4(k);
            p.m4*force.A_Y_G4(k);
            I_O4*acc.theta4_doubledot(k) ...
            ];

        x = A\B;

        force.F12x(k) = x(1);
        force.F12y(k) = x(2);
        force.F32x(k) = x(3);
        force.F32y(k) = x(4);
        force.F34x(k) = x(5);
        force.F34y(k) = x(6);
        force.F41x(k) = x(7);
        force.F41y(k) = x(8);
        force.M12(k)  = x(9);
    end

    % Convert forces to N and moments to both N*cm and N*m
    force.F12x_N = force.F12x./p.N_to_gcm;
    force.F12y_N = force.F12y./p.N_to_gcm;
    force.F32x_N = force.F32x./p.N_to_gcm;
    force.F32y_N = force.F32y./p.N_to_gcm;
    force.F34x_N = force.F34x./p.N_to_gcm;
    force.F34y_N = force.F34y./p.N_to_gcm;
    force.F41x_N = force.F41x./p.N_to_gcm;
    force.F41y_N = force.F41y./p.N_to_gcm;

    force.M12_Ncm = force.M12./p.N_to_gcm;
    force.M12_Nm  = force.M12_Ncm.*p.Ncm_to_Nm;

    % Shaking force  [N] and shaking moment [N*cm & N*m]
    force.Fs_x_N = -force.F12x_N + force.F41x_N;
    force.Fs_y_N = -force.F12y_N + force.F41y_N;
    force.Fs_mag_N = sqrt(force.Fs_x_N.^2 + force.Fs_y_N.^2);
    force.Fs_dir_rad = atan2(force.Fs_y_N, force.Fs_x_N);

    force.Ms_Ncm = force.M12_Ncm + force.F41y_N.*p.r1;
    force.Ms_Nm  = force.Ms_Ncm.*p.Ncm_to_Nm;

    % Turning-pair force magnitudes [N]
    force.F12_mag_N = sqrt(force.F12x_N.^2 + force.F12y_N.^2);
    force.F23_mag_N = sqrt(force.F32x_N.^2 + force.F32y_N.^2);
    force.F34_mag_N = sqrt(force.F34x_N.^2 + force.F34y_N.^2);
    force.F14_mag_N = sqrt(force.F41x_N.^2 + force.F41y_N.^2);
end

function plotPartB(p, pos, vel, acc, cp)
    %%  B.5 plots
    figure('Name', 'B.5(a) Angular Displacements');
    plot(p.theta2_deg, pos.theta3_deg);
    hold on;
    plot(p.theta2_deg, pos.theta4_deg);
    markPoint60(p, pos.theta3_deg);
    markPoint60(p, pos.theta4_deg);
    title('Angular Displacements');
    xlabel('\theta_2 [deg]');
    ylabel('\theta_3 and \theta_4 [deg]');
    legend({'\theta_3','\theta_4'}, 'Location', 'southeast');
    grid on;
    hold off;

    figure('Name', 'B.5(b) Angular Velocities');
    plot(p.theta2_deg, vel.theta3_dot);
    hold on;
    plot(p.theta2_deg, vel.theta4_dot);
    markPoint60(p, vel.theta3_dot);
    markPoint60(p, vel.theta4_dot);
    title('Angular Velocities');
    xlabel('\theta_2 [deg]');
    ylabel('\omega_3 and \omega_4 [rad/s]');
    legend({'\omega_3','\omega_4'}, 'Location', 'southeast');
    grid on;
    hold off;

    figure('Name', 'B.5(c) Angular Accelerations');
    plot(p.theta2_deg, acc.theta3_doubledot);
    hold on;
    plot(p.theta2_deg, acc.theta4_doubledot);
    markPoint60(p, acc.theta3_doubledot);
    markPoint60(p, acc.theta4_doubledot);
    title('Angular Accelerations');
    xlabel('\theta_2 [deg]');
    ylabel('\alpha_3 and \alpha_4 [rad/s^2]');
    legend({'\alpha_3','\alpha_4'}, 'Location', 'southeast');
    grid on;
    hold off;

    figure('Name', 'B.5(d) Cartesian Position of P');
    plot(p.theta2_deg, cp.xp);
    hold on;
    plot(p.theta2_deg, cp.yp);
    markPoint60(p, cp.xp);
    markPoint60(p, cp.yp);
    title('Cartesian Position Components of P');
    xlabel('\theta_2 [deg]');
    ylabel('x_P and y_P [cm]');
    legend({'x_P','y_P'}, 'Location', 'southeast');
    grid on;
    hold off;

    figure('Name', 'B.5(e) Velocity Magnitude of P');
    plot(p.theta2_deg, cp.vp_mag);
    hold on;
    markPoint60(p, cp.vp_mag);
    title('Velocity Magnitude of P');
    xlabel('\theta_2 [deg]');
    ylabel('|v_P| [cm/s]');
    legend({'|v_P|'}, 'Location', 'southeast');
    grid on;
    hold off;

    figure('Name', 'B.5(f) Acceleration Magnitude of P');
    plot(p.theta2_deg, cp.ap_mag);
    hold on;
    markPoint60(p, cp.ap_mag);
    title('Acceleration Magnitude of P');
    xlabel('\theta_2 [deg]');
    ylabel('|a_P| [cm/s^2]');
    legend({'|a_P|'}, 'Location', 'southeast');
    grid on;
    hold off;

    figure('Name', 'B.5(g) Coupler Curve');
    plot(cp.xp, cp.yp);
    title('Coupler Curve: y_P vs. x_P');
    xlabel('x_P [cm]');
    ylabel('y_P [cm]');
    grid on;
end

function plotPartC(p, force)
    %% Part C.3 plots: torque, shaking force, shaking moment, joint forces
    % These plots use the loaded case, Fext = 0.5 N

    figure('Name', 'C.3(i) Driving Torque');
    plot(p.theta2_deg, force.M12_Ncm);
    title('Driving Torque: M_{12} vs. \theta_2');
    xlabel('\theta_2 [deg]');
    ylabel('M_{12} [N*cm]');
    grid on;

    figure('Name', 'C.3(ii) Polar Shaking Force');
    polarplot(force.Fs_dir_rad, force.Fs_mag_N);
    title('Shaking Force F_S [N]');
    legend({'F_S [N]'});
    grid on;

    figure('Name', 'C.3(iii) Shaking Force Components');
    plot(p.theta2_deg, force.Fs_x_N, 'b', p.theta2_deg, force.Fs_y_N, 'r');
    title('Shaking Force Components vs. \theta_2');
    xlabel('\theta_2 [deg]');
    ylabel('Force Components [N]');
    legend({'(F_S)_x [N]', '(F_S)_y [N]'}, 'Location', 'southeast');
    grid on;

    figure('Name', 'C.3(iv) Shaking Moment');
    plot(p.theta2_deg, force.Ms_Ncm);
    title('Shaking Moment M_S vs. \theta_2');
    xlabel('\theta_2 [deg]');
    ylabel('M_S [N*cm]');
    grid on;

    figure('Name', 'C.3(v) Turning Pair Forces');
    plot(p.theta2_deg, force.F12_mag_N, 'b', ...
         p.theta2_deg, force.F23_mag_N, 'g', ...
         p.theta2_deg, force.F34_mag_N, 'r', ...
         p.theta2_deg, force.F14_mag_N, 'y');
    title('Turning Pair Force Magnitudes vs. \theta_2');
    xlabel('\theta_2 [deg]');
    ylabel('Turning Pair Force Magnitude [N]');
    legend({'|F_{12}| [N]', '|F_{23}| [N]', '|F_{34}| [N]', ...
        '|F_{14}| [N]'}, 'Location', 'northeast');
    grid on;
end

function C4 = summarizeC4(p, forceLoaded, forceNoLoad)
    %% Part C.4 numerical summary values
    % (a) Get max moment and it's index
    [C4.maxTorque_Nm, idxTorque] = max(abs(forceLoaded.M12_Nm));

    % Convert for multiple units
    C4.maxTorque_Ncm = C4.maxTorque_Nm/p.Ncm_to_Nm;
    C4.thetaTorque_deg = p.theta2_deg(idxTorque);

    % (b) Get max shaking force and it's index
    [C4.maxFs_N, idxFs] = max(forceLoaded.Fs_mag_N);
    C4.thetaFs_deg = p.theta2_deg(idxFs);

    %(c) Arrange all forces within one and find the max of it and its index
    pairForces = [forceLoaded.F12_mag_N(:), forceLoaded.F23_mag_N(:), ...
                  forceLoaded.F34_mag_N(:), forceLoaded.F14_mag_N(:)];
    pairNames = {'F12', 'F23', 'F34', 'F14'};
    % Find and set the max pair force and column
    [C4.maxPairForce_N, linearIdx] = max(pairForces(:));
    [idxPair, pairCol] = ind2sub(size(pairForces), linearIdx);

    % Get the name of the column and degree from row
    C4.maxPairName = pairNames{pairCol};
    C4.thetaPair_deg = p.theta2_deg(idxPair);

    % (d): ran once with both 0.5N and 0.0N applied (torque vs
    % torquenoload)
    % Get max torque no load (as (a) already did torque with force)
    [C4.maxTorqueNoLoad_Nm, idxNoLoadTorque] = ...
        max(abs(forceNoLoad.M12_Nm));
    C4.maxTorqueNoLoad_Ncm = C4.maxTorqueNoLoad_Nm/p.Ncm_to_Nm;

    % find theta at found index above
    C4.thetaNoLoadTorque_deg = p.theta2_deg(idxNoLoadTorque);

    C4.torqueChange_Nm = C4.maxTorque_Nm - C4.maxTorqueNoLoad_Nm;
    C4.torqueChange_Ncm = C4.maxTorque_Ncm - C4.maxTorqueNoLoad_Ncm;
    C4.percentChange = 100*C4.torqueChange_Nm/C4.maxTorqueNoLoad_Nm;
end

function printC4Summary(C4)
    %% Print Part C.4 answers to the Command Window

    fprintf('\n================ PART C.4  ====================\n');
    fprintf(['(a) Max |M12| = %.3f N*m = %.3f N*cm at theta2 = ' ...
        '%.0f deg.\n'], ...
        C4.maxTorque_Nm, C4.maxTorque_Ncm, C4.thetaTorque_deg);

    fprintf(['(b) Max shaking force |FS| = %.3f N at theta2 = ' ...
        '%.0f deg.\n'], ...
        C4.maxFs_N, C4.thetaFs_deg);

    fprintf(['(c) ' ...
        'Largest turning-pair force = |%s| = ' ...
        '%.4f N at theta2 = %.0f deg.\n'], ...
        C4.maxPairName, C4.maxPairForce_N, C4.thetaPair_deg);

    fprintf(['(d) With Fext = 0 N, max |M12| = %.5f N*m = ' ...
        '%.4f N*cm at theta2 = %.0f deg.\n'], ...
        C4.maxTorqueNoLoad_Nm, C4.maxTorqueNoLoad_Ncm, ...
        C4.thetaNoLoadTorque_deg);
    fprintf(['Motor Torque Change from no-load to Fext = ' ...
        '0.5 N: %.5f N*m = %.4f N*cm, or %.2f%%.\n'], ...
        C4.torqueChange_Nm, C4.torqueChange_Ncm, C4.percentChange);
    fprintf('')
    fprintf('=======================================================\n\n');
end

function markPoint60(p, yData)
    %% Adds a red marker and coords label to the active axs at theta2=60
    x_deg_60 = p.theta2_deg(p.indexAt60);
    y_value_60 = yData(p.indexAt60);
    plot(x_deg_60, y_value_60, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    text(x_deg_60, y_value_60, sprintf('  (%.0f, %.3f)', x_deg_60, y_value_60));
end
