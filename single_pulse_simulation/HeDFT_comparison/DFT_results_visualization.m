%% HeDFT Visualization: 9 Angström
clear; close all;

% Paths
base9 = 'T:\github synchronized\Iodine_Helium_Simulation\single_pulse_simulation\HeDFT_comparison\9Angström';

% Load Data
[t9_v1, v9_1] = importfile_v2(fullfile(base9, 'data_vabs.csv')); % Adjust filename if 'v1.csv' exists
[t9_v2, v9_2] = importfile_v2(fullfile(base9, 'data_vabs2.csv'));
[t9_R, R9]    = importfile_R1_R2(fullfile(base9, 'R1-R2.csv'));

figure('Color', 'w', 'Name', 'HeDFT 9A Reference');
tlo9 = tiledlayout(2, 2, 'TileSpacing', 'compact');


% Plot 1: Velocity Atom 1
nexttile;
plot(t9_v1, abs(v9_1), 'LineWidth', 1.5, 'Color', '#0072BD');
title('Velocity: Iodine 1 (9Å)'); ylabel('v / Å/ps'); grid on;
xlim([0, 12]); % Set a standard window for comparison

% Plot 2: Velocity Atom 2
nexttile;
plot(t9_v2, abs(v9_2), 'LineWidth', 1.5, 'Color', '#D95319');
title('Velocity: Iodine 2 (9Å)'); grid on;
xlim([0, 12]); % Set a standard window for comparison

% Plot 3: Internuclear Distance
nexttile([1 2]); % Span across bottom row
plot(t9_R, R9, 'k', 'LineWidth', 2);
title('Internuclear Distance R (9Å)'); xlabel('t / ps'); ylabel('R / Å'); grid on;

 




%% HeDFT Visualization: 18 Angström
% Paths
martiPath = 'T:\Cloud\MATLAB iodine\MartiPiNotes\impurities_dynamics\impurities_dynamics\';

% Load Data
V1_18 = importfile_marti(fullfile(martiPath, 'vimp.I_1'));
V2_18 = importfile_marti(fullfile(martiPath, 'vimp.I_2'));
R1_18 = importfile_marti(fullfile(martiPath, 'rimp.I_1'));
R2_18 = importfile_marti(fullfile(martiPath, 'rimp.I_2'));

% Calculations
v18_1 = sqrt(sum(V1_18.vec.^2, 2));
v18_2 = sqrt(sum(V2_18.vec.^2, 2));
R18   = sqrt(sum((R1_18.vec - R2_18.vec).^2, 2));

figure('Color', 'w', 'Name', 'HeDFT 18A Reference');
tlo18 = tiledlayout(2, 2, 'TileSpacing', 'compact');

% Plot 1: Velocity Atom 1
nexttile;
plot(V1_18.t, v18_1, 'LineWidth', 1.5, 'Color', '#0072BD');
title('Velocity: Iodine 1 (18Å)'); ylabel('v / Å/ps'); grid on;

% Plot 2: Velocity Atom 2
nexttile;
plot(V2_18.t, v18_2, 'LineWidth', 1.5, 'Color', '#D95319');
title('Velocity: Iodine 2 (18Å)'); grid on;

% Plot 3: Internuclear Distance
nexttile([1 2]);
plot(R1_18.t, R18, 'k', 'LineWidth', 2);
title('Internuclear Distance R (18Å)'); xlabel('t / ps'); ylabel('R / Å'); grid on;