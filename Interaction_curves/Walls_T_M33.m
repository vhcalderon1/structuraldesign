% ======================================================================================
% T-Wall Interaction Diagram Analysis and Visualization
% ======================================================================================
% Author: Victor Calderón (July 2019) 
% Updated: Jefferson De la Cuba (February 2025)
% ------------------------------------------------------------------------------
% This code performs moment-curvature analysis and generates interaction diagrams for reinforced 
% concrete T-walls according to ACI-318 specifications. The implementation calculates both nominal 
% and design strength envelopes while handling path configurations for external datasets and outputs.
%
% Inputs:
%   - Material properties: fpc, fy, Es (kg/cm² units)
%   - Geometric properties: b_1, b_2, h_1, h_2 (cm units)
%   - Reinforcement configuration: As (cm²), eid (cm spacing)
%   - External dataset: T_M33.txt (kN & kN-m units in 2 columns)
%
% Outputs:
%   - 3 Interaction diagrams (PNG files)
%   - Console outputs of maximum axial capacity
%   - Generated figures saved to ../outputs/ directory
%

clear;
close all;
clc;

% Create outputs directory if non-existent (using relative path)
output_dir = fullfile('..', 'outputs');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

%% Interaction Diagram - T Wall
% Material Properties
fpc = 210;         % Concrete compressive strength (kg/cm²)
beta_c = 0.85;     % Concrete stress block factor
fy = 4200;         % Steel yield strength (kg/cm²)
Es = 2.0000e+06;   % Steel modulus of elasticity (kg/cm²)
ey = 0.0020;       % Steel yield strain (unitless)
ecu = 0.003;       % Ultimate concrete strain (per ACI-318)

% Geometric Properties (cm)
b_1 = 675;         % Flange width (cm)
b_2 = 10;          % Web width (cm)
h_1 = 10;          % Flange thickness (cm)
h_2 = 192;         % Total height (cm); note: web height = h_2 - h_1
Ag = b_1 * h_1 + b_2 * (h_2 - h_1);  % Gross cross-sectional area (cm²)

%% Reinforcement Distribution
n = 43;  % Number of reinforcement layers (from top to bottom)

% Steel areas for each layer (cm²)
As = [repmat(1.13,12,1); repmat(0.5,9,1); 5; ...
      repmat(0.5,9,1); repmat(1.13,12,1)];

% Spacing between bars for each layer (cm) (from top to bottom)
eid = [3; repmat(8,6,1); 12; repmat(20,13,1); 14; 14; repmat(20,13,1); 12; ...
       repmat(8,6,1)];

% Calculate cumulative distances of reinforcement layers from the top (cm)
ei = [];
for i = 1:length(eid)-1
    ei(1) = eid(1);
    ei = [ei, ei(i) + eid(i+1)];
end
ei = ei';

n_c = 1;  % Step increment for neutral axis depth

%% Compression Analysis

% Calculate the centroid of the gross section (from the top)
y_c_1 = b_1 / 2;
y_c_2 = b_1 - y_c_1;  % (Not used later, but computed for completeness)

% Define a range of neutral axis depth values (c) to evaluate (cm)
c = 0.00001:n_c:b_1 + 0.3 * b_1;

% Initialize arrays for storing strain distributions and reversed reinforcement data
e_i_a = [];  % To store strain distributions for each neutral axis depth
As_2 = [];
ei_2 = [];

% Reverse the reinforcement order for the opposite side analysis
for i_3 = 0:length(As)-1
    As_2 = [As_2; As(n - i_3)];
    ei_2 = [ei_2; b_1 - ei(n - i_3)];
end

flag = 0;  % Flag to control iterations

% Maximum axial compression load with stirrups (kN)
% Convert kg to kN: multiply by 0.00980665
P_n_max = 0.8 * (0.85 * fpc * (Ag - sum(As)) + fy * sum(As)) * 0.00980665;

% Limits for the compression block based on geometry (cm)
lim_1 = (b_1 - b_2) / 2;
lim_2 = (b_1 + b_2) / 2;

% Perform two iterations: first for one side and then for the opposite side
while flag == 0 || flag == 1
    % Initialize arrays for each evaluated value of neutral axis depth c
    a = zeros(length(c),1);         % Depth of equivalent concrete stress block (cm)
    P_n = zeros(length(c),1);         % Nominal axial load (kN)
    P_real = zeros(length(c),1);      % Actual axial load (kN)
    M_n = zeros(length(c),1);         % Nominal moment (kN-m)
    M_real = zeros(length(c),1);      % Actual moment (kN-m)
    phi_Mn = zeros(length(c),1);      % Design moment (kN-m)
    phi_Pn = zeros(length(c),1);      % Design axial load (kN)
    
    % Loop over each neutral axis depth value
    for i = 1:length(c)
        e_i = zeros(length(As),1);    % Strains in reinforcement layers
        f_s = zeros(length(As),1);    % Forces in reinforcement layers (kg)
        ms_i = zeros(length(As),1);   % Moment contributions from reinforcement layers (kg-cm)
        
        % Compute the depth of the equivalent concrete stress block
        a(i) = c(i) * beta_c;
        
        % Loop over each reinforcement layer to calculate strain and force
        for i_1 = 1:length(As)
            % Calculate strain in the reinforcement layer
            e_i(i_1) = ecu * (ei(i_1) - c(i)) / c(i);
            % (Negative strain = compression; Positive = tension)
            
            % Calculate force in the reinforcement layer (kg)
            if abs(e_i(i_1)) > ey
                if e_i(i_1) > 0
                    f_s(i_1) = fy * As(i_1);
                else
                    f_s(i_1) = -fy * As(i_1);
                end
            else
                f_s(i_1) = Es * e_i(i_1) * As(i_1);
            end
            
            % Calculate concrete compression force (kg)
            if a(i) < lim_1
                P_con = 0.85 * fpc * a(i) * h_1;
            elseif lim_1 <= a(i) && a(i) < lim_2
                P_con_1 = 0.85 * fpc * h_1 * (b_1 - b_2) / 2;
                P_con_2 = 0.85 * fpc * (a(i) - (b_1 - b_2) / 2) * h_2;
                P_con = P_con_1 + P_con_2;
            elseif a(i) >= lim_2
                P_con_1 = 0.85 * fpc * b_2 * (h_2 - h_1);
                P_con_2 = 0.85 * fpc * h_1 * a(i);
                P_con = P_con_1 + P_con_2;
            end
            
            % Calculate the moment contribution of the reinforcement layer (kg-cm)
            ms_i(i_1) = f_s(i_1) * (y_c_1 - ei(i_1));
        end
        
        % Store the strain distribution for this neutral axis depth
        e_i_a = [e_i_a; e_i'];
        
        % Calculate the nominal moment (M_n) (kN-m)
        % Convert kg-cm to kN-m: multiply by 0.0000980665
        if a(i) < lim_1
            M_real(i) = (P_con * (y_c_1 - a(i)/2) - sum(ms_i)) * 0.0000980665;
        elseif lim_1 <= a(i) && a(i) < lim_2
            M_real(i) = (P_con_1 * (y_c_1 - (b_1 - b_2)/4) + ...
                         P_con_2 * (y_c_1 - (a(i) + (b_1 - b_2)/2)/2) - ...
                         sum(ms_i)) * 0.0000980665;
        elseif a(i) >= lim_2
            M_real(i) = (P_con_2 * (y_c_1 - a(i)/2) - sum(ms_i)) * 0.0000980665;
        end
        
        % Calculate the actual axial load (P_real) (kN)
        % Convert kg to kN: multiply by 0.00980665
        P_real(i) = (P_con - sum(f_s)) * 0.00980665;
        % (Positive P_real indicates compression; negative indicates tension)
        
        % Ensure the axial load does not exceed the maximum allowable compression load
        if P_n_max >= P_real(i)
            P_n(i) = P_real(i);
            M_n(i) = M_real(i);
        else
            M_n(i) = 0;
            P_n(i) = P_n_max;
        end
    end
    
    % Additional checks based on beam properties
    % Convert kg to kN: multiply by 0.00980665
    P_viga = 0.1 * fpc * Ag * 0.00980665;  % (kN)
    M_b = max(M_n);
    p_1 = find(M_n == M_b);          % Index of maximum moment
    P_b = P_n(p_1);
    
    P_tran_1 = P_viga / 0.7;  % For walls with stirrups
    P_tran_2 = P_b / 0.7;
    P_tran = min(P_tran_1, P_tran_2);
    
    % Calculate design strengths (phi factors) for moment and axial load
    for i_5 = 1:length(c)
        phi = 0.9 - 0.2 * P_n(i_5) / P_tran;
        if phi <= 0.7
            phi = 0.7;
        elseif phi >= 0.9
            phi = 0.9;
        end
        phi_Mn(i_5) = M_n(i_5) * phi;
        phi_Pn(i_5) = P_n(i_5) * phi;
    end
    
    % Store data for negative and positive sides of the interaction diagram
    flag = flag + 1;
    if flag == 1  % First iteration: negative side
        M_n_neg = -M_n;
        P_n_neg = P_n;
        e_i_a_neg = e_i_a;
        
        phi_Mn_neg = -phi_Mn;
        phi_Pn_neg = phi_Pn;
    else  % Second iteration: positive side
        M_n_pos = M_n;
        P_n_pos = P_n;
        e_i_a_pos = e_i_a;
        
        phi_Mn_pos = phi_Mn;
        phi_Pn_pos = phi_Pn;
    end
end

% Reverse the positive side data to close the interaction diagram
M_n_pos_graf = [];
P_n_pos_graf = [];
phi_M_n_pos_graf = [];
phi_P_n_pos_graf = [];
n_3 = length(M_n_pos);

for i_4 = 0:n_3-1
    M_n_pos_graf = [M_n_pos_graf; M_n_pos(n_3-i_4)];
    P_n_pos_graf = [P_n_pos_graf; P_n_pos(n_3-i_4)];
    phi_M_n_pos_graf = [phi_M_n_pos_graf; phi_Mn_pos(n_3-i_4)];
    phi_P_n_pos_graf = [phi_P_n_pos_graf; phi_Pn_pos(n_3-i_4)];
end

% Combine negative and positive side data for a complete interaction diagram
M_n_graf = [M_n_neg; M_n_pos_graf];
P_n_graf = [P_n_neg; P_n_pos_graf];

phi_M_n_graf = [phi_Mn_neg; phi_M_n_pos_graf];
phi_P_n_graf = [phi_Pn_neg; phi_P_n_pos_graf];

%% Plotting Section

% Configure output path using path traversal
save_path = @(fname) fullfile('..', 'outputs', fname);

% Plot 1: Nominal Interaction Diagram
figure;
plot(M_n_pos, P_n_pos, 'b-', 'LineWidth', 2);
hold on; grid on;
plot(M_n_neg, P_n_neg, 'r-', 'LineWidth', 2);
xlabel('Moment (kN-m)'); 
ylabel('Axial Load (kN)');
title('Nominal Interaction Diagram');
legend('Positive Side', 'Negative Side', 'Location', 'best');
print('-dpng', '-r300');
close(gcf);

% Plot 2: Complete Interaction Diagrams
figure;
plot(M_n_graf, P_n_graf, 'b-', 'LineWidth', 2);
hold on; grid on;
plot(phi_M_n_graf, phi_P_n_graf, 'k--', 'LineWidth', 2);
xlabel('Moment (kN-m)'); 
ylabel('Axial Load (kN)');
title('Complete Interaction Diagrams');
legend('Nominal', 'Design', 'Location', 'best');
print('-dpng', '-r300');
close(gcf);

% Wall Load Demands

% Load load data from an external dataset
% (The file T_M33.txt must be located in the "datasets" folder one level up)
load('../datasets/T_M33.txt');
P_u_M33 = T_M33(:,1);  % Axial load demands (kN)
M_u_M33 = T_M33(:,2);  % Moment demands (kN-m)

sz = 15;  % Marker size for scatter plot

% Define axis limits for plotting (kN-m for moment, kN for axial load)
x_lim = 2.5e3;   % Limit for moment (kN-m)
y_lim_sup = 3e3; % Upper limit for axial load (kN)
y_lim_inf = 1e3; % Lower limit for axial load (kN)

% Plot 3: Final Diagram with Load Demands
figure;
plot(M_n_graf, P_n_graf, 'Color', [0 0.8 1], 'LineWidth', 2);
hold on; grid on;
plot(phi_M_n_graf, phi_P_n_graf, '--', 'Color', [0 0.8 1], 'LineWidth', 2);
scatter(M_u_M33, P_u_M33, 15, 'MarkerFaceColor', [0.7 0.7 0.7]);
xlabel('Moment (kN-m)'); 
ylabel('Axial Load (kN)');
title('Final Interaction Diagram with Load Demands');
legend('Nominal', 'Design', 'Load Data', 'Location', 'best');
print(save_path('Final_Interaction_Diagram'), '-dpng', '-r300');
close(gcf);


%% Calculation of Maximum Axial Load
% Convert kg to kN: multiply by 0.00980665
Pu_Nch = Ag * fpc * 0.3 * 0.00980665;  % Maximum axial load (kN)