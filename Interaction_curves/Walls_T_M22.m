% ===========================================================================
% T-Wall Interaction Diagram Generator: Structural Capacity Analysis
% ===========================================================================
% Author: Victor Calderón (July 2019)
% Updated: Jefferson De la Cuba (February 2025)
% --------------------------------------------------------------------------
% This code generates nominal and factored interaction diagrams for T-shaped
% structural walls per ACI 318 provisions. The program calculates axial
% load-moment capacity pairs considering concrete crushing and steel yielding
% limits, then compares results with applied load combinations.
%
% INPUTS:
% - Material: fpc = 210 kg/cm², fy = 4200 kg/cm²
% - Geometry: Wall dimensions [cm] (b1, b2, h1, h2)
% - Reinforcement: Layer quantities (n), areas (As), spacing (eid)
% - Dataset: T_P34_M22.txt (Applied loads: P_u_M22 [kN], M_u_M22 [kN-m])
%
% OUTPUTS:
% - Interaction diagrams (Nominal/Factored)
% - PNG figures saved to ../outputs/
% - Console output of maximum allowed axial load

clear
close all
clc

%% Configure Output Directory
[currPath, ~, ~] = fileparts(mfilename('fullpath'));
outputDir = fullfile(currPath, '..', 'outputs');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end
%% Interaction Diagram - T-Wall

% Input - Material Properties
fpc = 210;           % Concrete compressive strength (kg/cm²)
betha_c = 0.85;      % Concrete compression block factor (per ACI-318)
fy = 4200;           % Steel yield strength (kg/cm²)
Es = 2.0e6;          % Modulus of elasticity of steel (kg/cm²)
ey = 0.0020;         % Steel yield strain
ecu = 0.003;         % Ultimate concrete strain (ACI-318)

% Input - Geometric Properties (cm)
b_1 = 394;           % Flange width
b_2 = 10;            % Web width
h_1 = 12;            % Flange thickness (height)
h_2 = 412;           % Total height (cm) – (h₂ – h₁ gives the web height)
Ag = b_1 * h_1 + b_2 * (h_2 - h_1);  % Gross area (cm²)

%% Input - Reinforcement Distribution

n = 26;  % Number of reinforcement layers

% From top to bottom:
As = [10.5; repmat(0.5, 18, 1); repmat(1.13, 7, 1)];  % Reinforcement area (cm²)

eid = [6; 9; repmat(20, 17, 1); 6; repmat(8, 7, 1)];     % Spacing between bars (cm)

ei = [];  % Distance of each layer from the top (cm)
for i = 1:length(eid)-1
    % The first element is fixed; subsequent positions are cumulative
    ei(1) = eid(1);
    ei = [ei, ei(i) + eid(i+1)];
end
ei = ei';  % Convert to a column vector

n_c = 1;  % Increment for neutral axis depth

%% Compression Analysis

% Calculate the centroid of the gross section (from top)
y_c_1 = (1/2 * ((b_1 * h_1^2) + (h_2^2 - h_1^2) * b_2)) / (b_1 * h_1 + b_2 * (h_2 - h_1));
y_c_2 = h_2 - y_c_1;  % Distance from the bottom

% Define a range for the neutral axis depth 'c' (cm)
c = 0.00001:n_c:h_2;  % Modify this range as needed for finer resolution

e_i_a = [];  % Matrix to store strains in reinforcement for each neutral axis depth

% Create mirrored reinforcement arrays for the negative bending case
As_2 = [];
ei_2 = [];
for i_3 = 0:length(As)-1
    As_2 = [As_2; As(n - i_3)];
    ei_2 = [ei_2; h_2 - ei(n - i_3)];
end

flag = 0;

% Maximum compressive axial load capacity with ties (kNs)
P_n_max = 0.8 * (0.85 * fpc * (Ag - sum(As)) + fy * sum(As)) * 9.81 / 1000; % Corrected unit conversion

while flag == 0 || flag == 1
    a = zeros(length(c), 1);         % Depth of equivalent stress block (cm)
    P_n = zeros(length(c), 1);         % Nominal axial load (kN)
    P_real = zeros(length(c), 1);      % Actual axial force (kN)
    M_n = zeros(length(c), 1);         % Nominal moment (kN-m)
    M_real = zeros(length(c), 1);      % Actual moment (kN-m)
    phi_Mn = zeros(length(c), 1);      % Strength-reduced moment (kN-m)
    phi_Pn = zeros(length(c), 1);      % Strength-reduced axial load (kN)

    for i = 1:length(c)
        e_i = zeros(length(As), 1);    % Strain in each reinforcement layer
        f_s = zeros(length(As), 1);    % Force in each reinforcement layer (kg)
        ms_i = zeros(length(As), 1);   % Moment contribution from each layer (kg-cm)

        % Calculate the equivalent rectangular stress block depth
        a(i) = c(i) * betha_c;

        for i_1 = 1:length(As)
            % Calculate strain in reinforcement layer
            % (negative indicates compression, positive indicates tension)
            e_i(i_1) = ecu * (ei(i_1) - c(i)) / c(i);

            % Calculate force in reinforcement (kg)
            if abs(e_i(i_1)) > ey
                if e_i(i_1) > 0
                    f_s(i_1) = fy * As(i_1);
                else
                    f_s(i_1) = -fy * As(i_1);
                end
            else
                f_s(i_1) = Es * e_i(i_1) * As(i_1);
            end

            % Calculate concrete compressive force (kg)
            if a(i) < h_1
                P_con = 0.85 * fpc * a(i) * b_1;
            else
                P_con_ala = 0.85 * fpc * h_1 * b_1;
                P_con_alma = 0.85 * fpc * (a(i) - h_1) * b_2;
                P_con = P_con_ala + P_con_alma;
            end

            % Calculate moment contribution from this reinforcement layer (kg-cm)
            ms_i(i_1) = f_s(i_1) * (y_c_1 - ei(i_1));
        end

        % Store reinforcement strains for the current neutral axis depth
        e_i_a = [e_i_a; e_i'];

        % Calculate the nominal moment (kN-m)
        if a(i) < h_1
            M_real(i) = (P_con * (y_c_1 - a(i)/2) - sum(ms_i)) * 9.81 / 100000; % Corrected unit conversion
        else
            M_real(i) = (P_con_ala * (y_c_1 - h_1/2) + P_con_alma * (y_c_1 - (a(i) - h_1)/2) - sum(ms_i)) * 9.81 / 100000; % Corrected
        end

        % Calculate the total axial force (kN)
        P_real(i) = (P_con - sum(f_s)) * 9.81 / 1000;  % Corrected unit conversion
        
        if P_n_max >= P_real(i)
            P_n(i) = P_real(i);
            M_n(i) = M_real(i);
        else
            M_n(i) = 0;
            P_n(i) = P_n_max;
        end
    end

    % Calculate the axial load capacity for a beam (for ties) in kN
    P_viga = 0.1 * fpc * Ag * 9.81 / 1000; % Corrected unit conversion
    M_b = max(M_n);
    p_1 = find(M_n == M_b);  % Index for maximum moment
    P_b = P_n(p_1);

    P_tran_1 = P_viga / 0.7;  % Based on beam capacity with ties
    P_tran_2 = P_b / 0.7;

    P_tran = min(P_tran_1, P_tran_2);

    % Calculate strength reduction factors for moment and axial force
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

    % Store the results and switch reinforcement for the opposite bending axis
    flag = flag + 1;
    if flag == 0 || flag == 1
        % For negative bending (mirrored reinforcement layout)
        M_n_neg = -M_n;
        P_n_neg = P_n;
        e_i_a_neg = e_i_a;

        phi_Mn_neg = -phi_Mn;
        phi_Pn_neg = phi_Pn;

        % Switch reinforcement arrays for the negative bending case
        As = As_2;
        ei = ei_2;
        h_1_r = h_1;
        h_2_r = h_2;
        h_1 = h_2_r;
        h_2 = h_1_r;
        b_1_r = b_1;
        b_2_r = b_2;
        b_1 = b_2_r;
        b_2 = b_1_r;
        y_c_1 = y_c_2;
    else
        % For positive bending
        M_n_pos = M_n;
        P_n_pos = P_n;
        e_i_a_pos = e_i_a;

        phi_Mn_pos = phi_Mn;
        phi_Pn_pos = phi_Pn;
    end
end

% Reorder the positive bending diagram so that the interaction diagram closes
M_n_pos_graf = [];
P_n_pos_graf = [];
phi_M_n_pos_graf = [];
phi_P_n_pos_graf = [];
n_3 = length(M_n_pos);
for i_4 = 0:n_3-1
    M_n_pos_graf = [M_n_pos_graf; M_n_pos(n_3 - i_4)];
    P_n_pos_graf = [P_n_pos_graf; P_n_pos(n_3 - i_4)];
    phi_M_n_pos_graf = [phi_M_n_pos_graf; phi_Mn_pos(n_3 - i_4)];
    phi_P_n_pos_graf = [phi_P_n_pos_graf; phi_Pn_pos(n_3 - i_4)];
end

% Combine negative and positive bending results
M_n_graf = [M_n_neg; M_n_pos_graf];
P_n_graf = [P_n_neg; P_n_pos_graf];

phi_M_n_graf = [phi_Mn_neg; phi_M_n_pos_graf];
phi_P_n_graf = [phi_Pn_neg; phi_P_n_pos_graf];

%% Plotting the Interaction Diagrams

% Figure 1: Nominal Interaction Diagram (Positive and Negative Bending)
figure
plot(M_n_pos, P_n_pos, 'b-', 'LineWidth', 2)
hold on
plot(M_n_neg, P_n_neg, 'r-', 'LineWidth', 2)
grid on
xlabel('Moment (kN-m)', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('Axial Load (kN)', 'FontSize', 12, 'FontName', 'Times New Roman')
title('Nominal Interaction Diagram - T-Wall', 'FontSize', 14)
legend('Positive Bending', 'Negative Bending', 'Location', 'best')

% Figure 2: Combined Nominal and Strength-Reduced (Factored) Interaction Diagrams
figure
plot(M_n_graf, P_n_graf, 'Color', [0 0.8 1], 'LineWidth', 2)
hold on
plot(phi_M_n_graf, phi_P_n_graf, '--', 'Color', [0 0.8 1], 'LineWidth', 2)
grid on
xlabel('Moment (kN-m)', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('Axial Load (kN)', 'FontSize', 12, 'FontName', 'Times New Roman')
title('Combined Nominal and Factored Interaction Diagram - T-Wall', 'FontSize', 14)
legend('Nominal', 'Factored', 'Location', 'best')

%% T-Wall Applied Load Data

% Load applied load data from the datasets folder (relative path)
dataPath = fullfile('..', 'datasets', 'T_P34_M22.txt');
load(dataPath)  % This loads variable T_P34_M22
P_u_M22 = T_P34_M22(:,1);  % Applied axial loads (kN)
M_u_M22 = T_P34_M22(:,2);  % Applied moments (kN-m)
sz = 15;  % Size of scatter plot markers

x_lim = 1.5e4;   % Limit for moment axis (kN-m)
y_lim_sup = 1.4e4; % Upper limit for axial load (kN)
y_lim_inf = 2e3; % Lower limit for axial load (kN)

% Figure 3: Interaction Diagram with Applied Loads
figure
plot(M_n_graf, P_n_graf, 'Color', [0 0.8 1], 'LineWidth', 2)
hold on
plot(phi_M_n_graf, phi_P_n_graf, '--', 'Color', [0 0.8 1], 'LineWidth', 2)
scatter(M_u_M22, P_u_M22, 45, 'MarkerEdgeColor', [0.3 0.3 0.3],...
    'MarkerFaceColor', [0.7 0.7 0.7], 'LineWidth', 0.8)
plot([-x_lim, x_lim], [0, 0], 'k-', 'LineWidth', 1.5)
plot([0, 0], [-y_lim_inf, y_lim_sup], 'k-', 'LineWidth', 1.5)

% Axis formatting
xticks(-x_lim:5e3:x_lim)  % 2000 kN-m increments
yticks(-y_lim_inf:2e3:y_lim_sup)  % 2000 kN increments (includes negatives)
xlabel('Moment (kN-m)', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('Axial Load (kN)', 'FontSize', 14, 'FontName', 'Times New Roman')

% Disable scientific notation explicitly
ax = gca;
ax.XAxis.Exponent = 0;  % Turns off X-axis scientific notation
ax.YAxis.Exponent = 0;  % Turns off Y-axis scientific notation

% Adjust grid and frame
box off
grid on
set(gca, 'LineWidth', 1.3, 'FontSize', 12, 'FontName', 'Times New Roman')
axis([-x_lim x_lim -y_lim_inf y_lim_sup])
title('T-Wall Interaction Diagram with Applied Loads', 'FontSize', 14)
legend('Nominal', 'Factored', 'Applied Loads', 'Location', 'best')
%% Enhanced Figure Handling
figNames = {'InteractionDiagram_Nominal', 'InteractionDiagram_Factored',...
           'InteractionDiagram_AppliedLoads'};

% Save all open figures except the Nominal Interaction Diagram
figHandles = findobj('Type', 'figure');
for idx = 1:length(figHandles)
    currentFig = figHandles(idx);
    figNumber = currentFig.Number;
    % Skip saving the first figure (Nominal Interaction Diagram)
    if figNumber == 1 || figNumber == 2
        close(gcf)
        continue;
    end
    figure(currentFig)
    outputPath = fullfile(outputDir, sprintf('Figure%d_%s.png',...
                   figNumber, figNames{figNumber}));
    print(outputPath, '-dpng', '-r300')
    fprintf('Saved: %s\n', outputPath)
    close(gcf)  % Close figure after saving
end
%% Display Final Results
fprintf('\nAnalysis completed successfully\n');
fprintf('Output directory: %s\n', outputDir);

%% Calculate Maximum Allowed Axial Load

Pu_Nch = Ag * fpc * 0.3 * 9.81 / 1000;  % Corrected unit conversion
fprintf('Maximum allowed axial load: %.2f kN\n', Pu_Nch);