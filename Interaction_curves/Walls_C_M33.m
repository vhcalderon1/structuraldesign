% =========================================================================
% INTERACTION DIAGRAM GENERATOR FOR REINFORCED CONCRETE T-WALL SECTIONS
% =========================================================================
% 
% Author:  Victor Calderón (November 2019)
% Updated: Jefferson De la Cuba (February 2025) 
% ------------------------------------------------------------------------
% Academic/Professional Purpose: Generates nominal and design interaction 
% diagrams (P-M curves) for reinforced concrete T-walls per ACI 318. 
% Performs strain compatibility analysis and considers material nonlinearity.
%
% Inputs:
%   - Material Properties: 
%       • fpc = Concrete compressive strength [kg/cm²]
%       • fy = Steel yield strength [kg/cm²]
%       • Es = Steel elastic modulus [kg/cm²]
%   - Geometric Properties: 
%       • Cross-section dimensions (b_1, b_2, b_3, h_1, h_2, h_3) [cm]
%       • Reinforcement layout (As, eid arrays) [cm², cm]
%   - External Data: 
%       • C_M33.txt - Applied loads dataset [kN, kN-m] 
%         (Located in ../datasets directory)
% Outputs:
%   - PNG Figures: 
%       1. Nominal_Interaction_Diagram.png
%       2. Design_Interaction_Diagram.png 
%       3. Applied_Loads_Diagram.png
%   - Console: 
%       • Maximum axial capacity [kN]
%       • Validation messages for path handling


clear
close all
clc

%% Directory Configuration
[scriptPath, ~] = fileparts(mfilename('fullpath')); % Get current script path

% Configure output directory using path traversal
outputDir = fullfile(scriptPath, '..', 'outputs');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

%% Interaction Diagram - T-Wall

% Element Properties
fpc = 210;           % Concrete compressive strength (kg/cm²)
beta_c = 0.85;       % Compression block factor for concrete
fy = 4200;           % Steel yield strength (kg/cm²)
Es = 2.0e+06;        % Modulus of elasticity of steel (kg/cm²)
ey = 0.0020;         % Steel yield strain (unit strain)
ecu = 0.003;         % Ultimate concrete compressive strain (per ACI-318)

% Geometric Properties (cm)
% Element 1: Top rectangle (smaller dimension)
% Element 2: Web
% Element 3: Bottom rectangle (larger dimension)
b_1 = 130;           % Width of Element 1 (cm)
b_2 = 10;            % Width of the web (cm)
b_3 = 172;           % Width of Element 3 (cm)
h_1 = 10;            % Height of Element 1 (cm)
h_2 = 285;           % Height of the web (cm)
h_3 = 10;            % Height of Element 3 (cm)
Ag = b_1*h_1 + b_2*h_2 + b_3*h_3;  % Gross area (cm²)

%% Reinforcement Distribution

n = 19;  % Number of reinforcement layers

% From top to bottom:
As = [9.41; 1.13; repmat(0.5,15,1); 1.13; 10.41];   % Reinforcement area per layer (cm²)
eid = [5; 5; 3; repmat(20,14,1); 2; 5];                % Spacing between bars (cm)
ei = [];  % Vertical positions of reinforcement (cm)

% Calculate cumulative vertical distances for reinforcement layers
for i = 1:length(eid)-1
    ei(1) = eid(1);
    ei = [ei, ei(i) + eid(i+1)];
end
ei = ei';

n_c = 1; % Increment for neutral axis depth

%% Compression Analysis

% Calculate the centroid of the gross section (from top to bottom)
y_c_1 = ((b_1*h_1^2)/2 + b_2*h_2*(h_1 + h_2/2) + b_3*h_3*(h_1 + h_2 + h_3/2)) / ...
        (b_1*h_1 + b_2*h_2 + b_3*h_3);
y_c_2 = (h_1 + h_2 + h_3) - y_c_1;

% Define the range for the neutral axis depth 'c' (cm)
c = 0.00001:n_c:(h_1+h_2+h_3) + 0.6*(h_1+h_2+h_3);  % Adjust this range as needed

e_i_a = [];  % Array to store strains in reinforcement layers

% Prepare arrays for the reversed reinforcement distribution (for the negative side)
As_2 = [];
ei_2 = [];
As_1 = As;  % Backup of the original reinforcement distribution

for i_3 = 0:length(As)-1
    As_2 = [As_2; As(n - i_3)];
    ei_2 = [ei_2; (h_1 + h_2 + h_3) - ei(n - i_3)];
end

flag = 0;

% Maximum compression load with stirrups (kN)
force_conversion = 101.97; % Conversion factor from kgf to kN
moment_conversion = 10197; % Conversion factor from kgf·cm to kN·m

P_n_max = 0.8 * (0.85 * fpc * (Ag - sum(As)) + fy * sum(As)) / force_conversion;  % (kN)
lim_1 = h_1;
lim_2 = h_1 + h_2;

while flag == 0 || flag == 1
    % Initialize arrays for each neutral axis depth iteration
    a = zeros(length(c),1);
    P_n = zeros(length(c),1);
    P_real = zeros(length(c),1);
    M_n = zeros(length(c),1);
    M_real = zeros(length(c),1);
    phi_Mn = zeros(length(c),1);
    phi_Pn = zeros(length(c),1);

    for i = 1:length(c)
        e_i = zeros(length(As),1);
        f_s = zeros(length(As),1);
        ms_i = zeros(length(As),1);

        a(i) = c(i) * beta_c;

        for i_1 = 1:length(As)
            % Calculate strain in the reinforcement bars
            e_i(i_1) = ecu * (ei(i_1) - c(i)) / c(i);
            % Note: Negative indicates compression; positive indicates tension

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

            % Calculate the concrete compression force (kg)
            if a(i) < lim_1
                P_con = 0.85 * fpc * a(i) * b_1;
            elseif a(i) >= lim_1 && a(i) < lim_2
                P_con_1 = 0.85 * fpc * h_1 * b_1;
                P_con_2 = 0.85 * fpc * (a(i) - h_1) * b_2;
                P_con = P_con_1 + P_con_2;
            elseif a(i) >= lim_2
                P_con_1 = 0.85 * fpc * b_1 * h_1;
                P_con_2 = 0.85 * fpc * h_2 * b_2;
                P_con_3 = 0.85 * fpc * b_3 * (a(i) - (h_1 + h_2));
                P_con = P_con_1 + P_con_2 + P_con_3;
            end

            % Calculate the moment contribution from the reinforcement layer (kg-cm)
            ms_i(i_1) = f_s(i_1) * (y_c_1 - ei(i_1));
        end

        % Store the strains for this neutral axis depth
        e_i_a = [e_i_a; e_i'];

        % Calculate the nominal moment (M_n) in kN-m
        if a(i) < lim_1
            M_real(i) = (P_con * (y_c_1 - a(i)/2) - sum(ms_i)) / moment_conversion;
        elseif a(i) >= lim_1 && a(i) < lim_2
            M_real(i) = (P_con_1 * (y_c_1 - h_1/2) + P_con_2 * (y_c_1 - (a(i) + h_1)/2) - sum(ms_i)) / moment_conversion;
        elseif a(i) >= lim_2
            M_real(i) = (P_con_1 * (y_c_1 - h_1/2) + P_con_2 * (y_c_1 - (h_1 + h_2/2)) + ...
                         P_con_3 * (y_c_1 - (a(i) + h_1 + h_2)/2) - sum(ms_i)) / moment_conversion;
        end

        % Calculate the total axial force (P_n) in kN 
        % (compression is positive; tension is negative)
        P_real(i) = (P_con - sum(f_s)) / force_conversion;
        
        if P_n_max >= P_real(i)
            P_n(i) = P_real(i);
            M_n(i) = M_real(i);
        else
            M_n(i) = 0;
            P_n(i) = P_n_max;
        end
    end
    
    % Additional calculations for beam (viga) effects
    P_viga = 0.1 * fpc * Ag / force_conversion;  % (kN)
    M_b = max(M_n);
    p_1 = find(M_n == M_b);  % Index of maximum moment
    P_b = P_n(p_1);
    
    P_tran_1 = P_viga / 0.7;  % For walls with stirrups
    P_tran_2 = P_b / 0.7;
    
    P_tran = min(P_tran_1, P_tran_2);
    
    % Calculate phi factors and the design strengths
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
    
    % Store data from the positive and negative side iterations
    flag = flag + 1;
    
    if flag == 1  % First iteration: positive side
        M_n_pos = M_n;
        P_n_pos = P_n;
        e_i_a_pos = e_i_a;
        
        phi_Mn_pos = phi_Mn;
        phi_Pn_pos = phi_Pn;
        
        % Prepare for the negative side by reversing the reinforcement distribution and geometry
        As = As_2;
        ei = ei_2;
        h_1_r = h_1;
        h_2_r = h_2;
        h_3_r = h_3;
        h_1 = h_3_r;
        h_2 = h_2_r;
        h_3 = h_1_r;
        b_1_r = b_1;
        b_2_r = b_2;
        b_3_r = b_3;
        b_1 = b_3_r;
        b_2 = b_2_r;
        b_3 = b_1_r;
        y_c_1 = y_c_2;
        
    else  % Second iteration: negative side
        M_n_neg = -M_n;
        P_n_neg = P_n;
        e_i_a_neg = e_i_a;
        
        phi_Mn_neg = -phi_Mn;
        phi_Pn_neg = phi_Pn;
    end
end

% Reverse the positive side arrays for proper closure of the interaction diagram
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

% Concatenate the negative and positive parts of the diagram
M_n_graf = [M_n_neg; M_n_pos_graf];
P_n_graf = [P_n_neg; P_n_pos_graf];

phi_M_n_graf = [phi_Mn_neg; phi_M_n_pos_graf];
phi_P_n_graf = [phi_Pn_neg; phi_P_n_pos_graf];

%% Plotting the Interaction Diagrams

% Figure 1: Nominal Interaction Diagrams for Positive and Negative Sides
figure
plot(M_n_pos, P_n_pos, 'b-', 'LineWidth', 2, 'DisplayName', 'Positive Side')
hold on
plot(M_n_neg, P_n_neg, 'r-', 'LineWidth', 2, 'DisplayName', 'Negative Side')
grid on
xlabel('Moment (kN-m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Axial Load (kN)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Nominal Interaction Diagrams (Positive and Negative)', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('Location','Best')
set(gca, 'LineWidth', 1.3, 'FontSize', 12, 'FontName', 'Times New Roman')

print(fullfile(outputDir, 'Nominal_Interaction_Diagram'), '-dpng', '-r300')
fprintf('Saved: Nominal_Interaction_Diagram.png\n');
close(gcf)  % Prevent memory buildup

% Figure 2: Combined Nominal and Design Interaction Diagrams
figure
plot(M_n_graf, P_n_graf, 'b-', 'LineWidth', 2, 'DisplayName', 'Nominal Interaction Diagram')
hold on
plot(phi_M_n_graf, phi_P_n_graf, 'k--', 'LineWidth', 2, 'DisplayName', 'Design Interaction Diagram')
grid on
xlabel('Moment (kN-m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Axial Load (kN)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Combined Nominal and Design Interaction Diagrams', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('Location','Best')
set(gca, 'LineWidth', 1.3, 'FontSize', 12, 'FontName', 'Times New Roman')

print(fullfile(outputDir, 'Combined_Nominal_and_Design_Interaction_Diagrams'), '-dpng', '-r300')
fprintf('Saved: Combined_Nominal_and_Design_Interaction_Diagrams.png\n');
close(gcf)  % Prevent memory buildup

%% Applied Loads on the Wall

% Define relative paths for the dataset and outputs folders
datasetDir = fullfile('..', 'datasets');
outputDir  = fullfile('..', 'outputs');

% Load applied load data (assumed to be in the "datasets" folder)
data = load(fullfile(datasetDir, 'C_M33.txt'));
P_u_M33 = data(:,1);  % Applied axial loads (kN)
M_u_M33 = data(:,2);  % Applied moments (kN-m)
sz = 15;              % Marker size for the scatter plot

x_lim = 2*10^3;       % x-axis limit for moment (kN-m)
y_lim_sup = 2*10^3;   % Upper y-axis limit for axial load (kN)
y_lim_inf = 0.5*10^3; % Lower y-axis limit for axial load (kN)

figure
p1 = plot(M_n_graf, P_n_graf, 'Color', [0 0.8 1], 'LineWidth', 2, 'DisplayName', 'Nominal Interaction Diagram');
hold on
p2 = plot(phi_M_n_graf, phi_P_n_graf, '--', 'Color', [0 0.8 1], 'LineWidth', 2, 'DisplayName', 'Design Interaction Diagram');
hold on
p3 = scatter(M_u_M33, P_u_M33, sz, 'MarkerEdgeColor', [0.7 0.7 0.7], ...
    'MarkerFaceColor', [0.7 0.7 0.7], 'LineWidth', 0.01, 'DisplayName', 'Applied Loads');
hold on
% Plot zero axial load and zero moment reference lines
plot([-x_lim x_lim], [0 0], 'k', 'LineWidth', 0.9)
plot([0 0], [-y_lim_inf y_lim_sup], 'k', 'LineWidth', 0.9)
grid on
xticks(-x_lim:10^3:x_lim);
yticks(-y_lim_inf:0.5*10^3:y_lim_sup);
xlabel('Moment (kN-m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Axial Load (kN)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Interaction Diagram with Applied Loads', 'FontSize', 14, 'FontName', 'Times New Roman');
ytickformat('%.1f')  % Format tick labels on the y-axis
ax = gca;
ax.YAxis.Exponent = 3; % Set exponent for y-axis tick labels
ax.XAxis.Exponent = 3; % Set exponent for x-axis tick labels
box off
set(gca, 'LineWidth', 1.3, 'FontSize', 12, 'FontName', 'Times New Roman')
axis([-x_lim x_lim -y_lim_inf y_lim_sup])
legend('show', 'Location', 'Best')

print(fullfile(outputDir, 'Interaction_Diagram_with_Applied_Loads'), '-dpng', '-r300')
fprintf('Saved: Interaction_Diagram_with_Applied_Loads.png\n');
close(gcf)  % Prevent memory buildup

%% Calculation of Maximum Axial Load
Pu_Nch = Ag * fpc * 0.3 / force_conversion;  % Maximum axial load (kN)
%% Final Outputs
fprintf('\nMaximum Axial Capacity: %.2f kN\n', Pu_Nch);
fprintf('All outputs saved to: %s\n', outputDir);
