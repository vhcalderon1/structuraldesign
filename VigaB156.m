%% BEAM INTERACTION DIAGRAM GENERATOR

% Description:
%   Generates nominal and design interaction diagrams for reinforced concrete 
%   beams under combined axial load and bending moment. Follows ACI 318 
%   provisions for strength reduction factors and strain compatibility.
%
% Inputs:
%   - Material: fpc, beta_c, fy, Es, ey, ecu
%   - Geometric: b1, b2, h1, h2, Ag
%   - Reinforcement: n, As, eid, ei
%   - Analysis: n_c (neutral axis increment)
%
% Outputs:
%   - Figures: Nominal and design interaction diagrams (kN·m vs kN)
%   - Data arrays: M_n, P_n, phi_Mn, phi_Pn for positive/negative bending
%
% Units:
%   - Inputs in kg/cm² (strength), cm (dimensions)
%   - Outputs in kN (forces) and kN·m (moments)
%
% Authors:
%   - Original: Victor Calderón (April 2020)
%   - Updated: Jefferson De la Cuba (February 2025)
% =========================================================================

clear;
close all;
clc;

% Configure output directory using path traversal
currentDir = pwd;                          % Get execution directory
outputDir = fullfile(currentDir, '..', 'outputs');  % Go up 1 level + outputs

% Create outputs directory if non-existent
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

%% Input - Material Properties

fpc = 210;          % Concrete compressive strength (kg/cm^2)
beta_c = 0.85;      % Compression block factor
fy = 4200;          % Steel yield strength (kg/cm^2)
Es = 2.0000e+06;    % Steel modulus of elasticity (kg/cm^2)
ey = 0.0020;        % Steel yield strain
ecu = 0.003;        % Ultimate concrete strain (per ACI-318)

%% Input - Geometric Properties (cm)

b1 = 0;             % Flange width
b2 = 50;            % Web width
h1 = 0;             % Flange height
h2 = 90;            % Total height (or overall height if flange is zero)
Ag = b1*h1 + b2*(h2 - h1);  % Gross area (cm^2)

%% Input - Reinforcement Distribution

n = 3;              % Number of reinforcement layers
As = [11.36; 10.2; 15.88];  % Steel area (cm^2)
eid = [6; 72; 6];   % Spacing between bars (cm)

ei = [];  % Distance from the top to each bar (cm)

% Calculate cumulative distances for the reinforcement layers
for i = 1:length(eid)-1
    ei(1) = eid(1);
    ei = [ei, ei(i) + eid(i+1)];
end

ei = ei';

n_c = 1;  % Step increment for neutral axis depth

%% Compression Analysis

% Calculate the centroid of the gross section (from top to bottom)
y_c1 = (1/2 * ((b1*h1^2) + (h2^2 - h1^2)*b2)) / (b1*h1 + b2*(h2 - h1));
y_c2 = h2 - y_c1;

% Define a vector of neutral axis depths to evaluate (cm)
c = 0.00001:n_c:h2;  % (Adjust as needed for resolution)

e_i_a = [];  % To store the strains at each neutral axis depth

% Prepare arrays for reversed reinforcement layers (for negative bending)
As_2 = [];
ei_2 = [];

As_1 = As;  % Original reinforcement distribution

for i_3 = 0:length(As)-1
    As_2 = [As_2; As(n - i_3)];
    ei_2 = [ei_2; h2 - ei(n - i_3)];
end

flag = 0;

% Maximum axial compression load (kN) - with stirrups 
P_n_max = 0.8 * (0.85 * fpc * (Ag - sum(As)) + fy * sum(As)) * 9.81 / 1000;

while flag == 0 || flag == 1
    % Initialize arrays for each iteration over neutral axis depths
    a       = zeros(length(c), 1);    % Depth of equivalent stress block (cm)
    P_n     = zeros(length(c), 1);    % Nominal axial load (kN)
    P_real  = zeros(length(c), 1);    % Actual axial load (kN)
    M_n     = zeros(length(c), 1);    % Nominal moment (kN.m)
    M_real  = zeros(length(c), 1);    % Actual moment (kN.m)
    phi_Mn  = zeros(length(c), 1);    % Reduced moment capacity (kN.m)
    phi_Pn  = zeros(length(c), 1);    % Reduced axial load capacity (kN)

    % Loop over each considered neutral axis depth
    for i = 1:length(c)
        e_i  = zeros(length(As), 1);   % Strain in each reinforcement bar
        f_s  = zeros(length(As), 1);   % Force in each reinforcement bar (kg)
        ms_i = zeros(length(As), 1);   % Moment contribution from each reinforcement layer (kg-cm)

        % Calculate depth of the equivalent concrete stress block
        a(i) = c(i) * beta_c;

        % Loop over each reinforcement layer
        for i_1 = 1:length(As)
            % Calculate strain in the reinforcement bar
            e_i(i_1) = ecu * (ei(i_1) - c(i)) / c(i);
            
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

            % Calculate concrete compressive force (kg)
            if a(i) < h1
                P_con = 0.85 * fpc * a(i) * b1;
            else
                P_con_flange = 0.85 * fpc * h1 * b1;
                P_con_web    = 0.85 * fpc * (a(i) - h1) * b2;
                P_con        = P_con_flange + P_con_web;
            end

            % Calculate the moment contribution from the reinforcement layer (kg-cm)
            ms_i(i_1) = f_s(i_1) * (y_c1 - ei(i_1));
        end

        % Store the strains for each neutral axis depth
        e_i_a = [e_i_a; e_i'];

        % Calculate the nominal moment M_n (kN.m) 
        if a(i) < h1
            M_real(i) = (P_con * (y_c1 - a(i)/2) - sum(ms_i)) * 9.81 / 100000;
        else
            M_real(i) = (P_con_flange * (y_c1 - h1/2) + P_con_web * (y_c1 - (a(i) - h1)/2) - sum(ms_i)) * 9.81 / 100000;
        end

        % Calculate the total axial force P_n (kN) 
        P_real(i) = (P_con - sum(f_s)) * 9.81 / 1000;
        
        if P_n_max >= P_real(i)
            P_n(i) = P_real(i);
            M_n(i) = M_real(i);
        else
            M_n(i) = 0;
            P_n(i) = P_n_max;
        end
    end
    
    % Calculate the shear capacity of the beam (kN) 
    P_beam = 0.1 * fpc * Ag * 9.81 / 1000;
    M_b = max(M_n);
    p_1 = find(M_n == M_b);  % Index corresponding to maximum moment
    P_b = P_n(p_1);
    
    % Calculate transformation axial load values (kN) for stirrups case
    P_tran_1 = P_beam / 0.7;
    P_tran_2 = P_b / 0.7;
    
    P_tran = min(P_tran_1, P_tran_2);
    
    % Calculate reduced capacities using strength reduction factors
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
    
    % Store the results and switch reinforcement configuration
    flag = flag + 1;
    
    if flag == 0 || flag == 1
        % Store negative bending side results
        M_n_neg = -M_n;
        P_n_neg = P_n;
        e_i_a_neg = e_i_a;
        
        phi_Mn_neg = -phi_Mn;
        phi_Pn_neg = phi_Pn;
        
        % Switch reinforcement configuration for negative bending
        As = As_2;
        ei = ei_2;
        h1_r = h1;
        h2_r = h2;
        h1 = h2_r;
        h2 = h1_r;
        b1_r = b1;
        b2_r = b2;
        b1 = b2_r;
        b2 = b1_r;
        y_c1 = y_c2;
    else
        % Store positive bending side results
        M_n_pos = M_n;
        P_n_pos = P_n;
        e_i_a_pos = e_i_a;
        
        phi_Mn_pos = phi_Mn;
        phi_Pn_pos = phi_Pn;
    end
end

% Reverse the order of the positive side arrays
M_n_pos_plot = [];
P_n_pos_plot = [];
phi_M_n_pos_plot = [];
phi_P_n_pos_plot = [];
n_3 = length(M_n_pos);

for i_4 = 0:n_3-1
    M_n_pos_plot = [M_n_pos_plot; M_n_pos(n_3 - i_4)];
    P_n_pos_plot = [P_n_pos_plot; P_n_pos(n_3 - i_4)];
    phi_M_n_pos_plot = [phi_M_n_pos_plot; phi_Mn_pos(n_3 - i_4)];
    phi_P_n_pos_plot = [phi_P_n_pos_plot; phi_Pn_pos(n_3 - i_4)];
end

% Combine negative and positive sides
M_n_plot = [M_n_neg; M_n_pos_plot];
P_n_plot = [P_n_neg; P_n_pos_plot];

phi_M_n_plot = [phi_Mn_neg; phi_M_n_pos_plot];
phi_P_n_plot = [phi_Pn_neg; phi_P_n_pos_plot];

%% Improved Plotting

% Plot Nominal Interaction Diagram
figure;
plot(M_n_pos, P_n_pos, 'b-', 'LineWidth', 2);
hold on;
plot(M_n_neg, P_n_neg, 'r-', 'LineWidth', 2);
xlabel('Nominal Moment, M_n (kN·m)');
ylabel('Nominal Axial Load, P_n (kN)');
title('Nominal Interaction Diagram');
legend('Positive Bending', 'Negative Bending', 'Location', 'Best');
grid on;

% Save the figure to outputs directory
nominalFigPath = fullfile(outputDir, 'Nominal_Interaction_Diagram.png');
saveas(gcf, nominalFigPath);
fprintf('Saved nominal interaction diagram to: %s\n', nominalFigPath);

% Close the figure
close(gcf);

% Plot Design Interaction Diagram
figure;
plot(M_n_plot, P_n_plot, 'b-', 'LineWidth', 2);
hold on;
plot(phi_M_n_plot, phi_P_n_plot, 'r--', 'LineWidth', 2);
xlabel('Moment, M (kN·m)');
ylabel('Axial Load, P (kN)');
title('Interaction Diagram with Strength Reduction');
legend('Nominal Strength', 'Design Strength (\phiM_n, \phiP_n)', 'Location', 'Best');
grid on;

% Save the figure to outputs directory
designFigPath = fullfile(outputDir, 'Design_Interaction_Diagram.png');
saveas(gcf, designFigPath);
fprintf('Saved design interaction diagram to: %s\n', designFigPath);

% Close the figure
close(gcf);