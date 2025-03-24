%% Interaction Diagram Analysis for Reinforced Concrete T-Beams
% 
% Author: Victor Calder�n (April 2020)
% Updated: Jefferson De la Cuba (February 2025)
% -------------------------------------------------------------------------
% This code generates an interaction diagram for reinforced concrete T-beams
% following ACI-318 standards. It calculates both nominal and strength-reduced
% (phi-factor) capacity curves considering material nonlinearities and
% combined axial-flexural loading conditions.
%
% Inputs:
%   Material Properties:
%   - fpc: Concrete compressive strength [kg/cm�]
%   - fy: Steel yield strength [kg/cm�]
%   - Es: Steel elastic modulus [kg/cm�]
%   - ecu: Ultimate concrete strain [-]
%   
%   Geometrical Properties:
%   - b_1/b_2: Flange/web widths [cm]
%   - h_1/h_2: Flange height/total height [cm]
%   
%   Reinforcement:
%   - n: Number of steel layers
%   - As: Steel areas array [cm�]
%   - eid: Reinforcement spacing array [cm]
%
% Outputs:
%   - Interaction diagrams (PNG format) saved in '/outputs' directory
%   - Console output of key capacity values (P_n_max, M_b)
% -------------------------------------------------------------------------

clear
close all
clc

% Configure output directory using path traversal
currentDir = pwd;                          % Get execution directory
outputDir = fullfile(currentDir, '..', 'outputs');  % Go up 1 level + outputs

% Create outputs directory if non-existent
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end
    
%% Interaction Diagram - T-beam
% Input - Material Properties
fpc = 210; % Concrete compressive strength (kg/cm�)
betha_c = 0.85; % Concrete compression block factor
fy = 4200; % Yield strength of steel (kg/cm�)
Es = 2.0000e+06; % Modulus of elasticity for steel (kg/cm�)
ey = 0.0020; % Steel strain at yield
ecu = 0.003; % Ultimate strain in concrete (according to ACI-318)

% Input - Geometrical Properties (cm)
b_1 = 0; % Width of the flange
b_2 = 50; % Width of the web
h_1 = 0; % Height of the flange
h_2 = 80; % Total height - Height of the web (h_2 - h_1)
Ag = b_1 * h_1 + b_2 * (h_2 - h_1); % Gross area of the section (cm�)

%% Input - Reinforcement Distribution
n = 3; % Number of steel layers
% n = 2; % Uncomment this for a 2-layer configuration

% Steel areas (cm�) from top to bottom
As = [15.88; 5.68; 11.36]; 
% As = [15.88; 11.36]; % Uncomment for different sections

eid = [6; 5; 63]; % Distance between reinforcement bars (cm)
% eid = [6; 68]; % Uncomment for different bar spacing

ei = []; % Empty array for distances

% Calculate the distances between reinforcement bars
for i = 1:length(eid) - 1
    ei(1) = eid(1);
    ei = [ei ei(i) + eid(i + 1)];
end
ei = ei';

n_c = 1; % Neutral axis advance

%% Compression

% Calculation of the centroid of the gross section
y_c_1 = (1 / 2 * ((b_1 * h_1^2) + (h_2^2 - h_1^2) * b_2)) / (b_1 * h_1 + b_2 * (h_2 - h_1));
y_c_2 = h_2 - y_c_1;

c = 0.00001 : n_c : h_2 + 0 * h_2; % You can adjust the number of points to observe
e_i_a = [];

As_2 = [];
ei_2 = [];

As_1 = As;

% Flip the reinforcement areas and distances for calculation
for i_3 = 0 : length(As) - 1
    As_2 = [As_2; As(n - i_3)];
    ei_2 = [ei_2; h_2 - ei(n - i_3)];
end

flag = 0;

% Maximum compressive load with stirrups (kNs)
P_n_max = 0.8 * (0.85 * fpc * (Ag - sum(As)) + fy * sum(As)) * 0.00980665;

% Main loop to calculate interaction diagram
while flag == 0 || flag == 1
    
    % Initialize arrays for results
    a = zeros(length(c), 1);
    P_n = zeros(length(c), 1);
    P_real = zeros(length(c), 1);
    M_n = zeros(length(c), 1);
    M_real = zeros(length(c), 1);
    phi_Mn = zeros(length(c), 1);
    phi_Pn = zeros(length(c), 1);
    
    for i = 1:length(c)
        e_i = zeros(length(As), 1);
        f_s = zeros(length(As), 1);
        ms_i = zeros(length(As), 1);

        a(i) = c(i) * betha_c;

        for i_1 = 1:length(As)

            % Calculate strain in the steel bars
            e_i(i_1) = ecu * (ei(i_1) - c(i)) / c(i);

            % Calculate the force in each steel layer (kg)
            if abs(e_i(i_1)) > ey
                if e_i(i_1) > 0
                    f_s(i_1) = fy * As(i_1);
                else
                    f_s(i_1) = -fy * As(i_1);
                end
            else
                f_s(i_1) = Es * e_i(i_1) * As(i_1);
            end

            % Calculate the concrete compressive force (kg)
            if a(i) < h_1
                P_con = 0.85 * fpc * a(i) * b_1;
            else
                P_con_ala = 0.85 * fpc * h_1 * b_1;
                P_con_alma = 0.85 * fpc * (a(i) - h_1) * b_2;
                P_con = P_con_ala + P_con_alma;
            end

            % Calculate moment generated by steel layers (kg�cm)
            ms_i(i_1) = f_s(i_1) * (y_c_1 - ei(i_1));
        end
        
        % Store strains for each neutral axis position
        e_i_a = [e_i_a; e_i'];

        % Calculate the nominal moment - Mn (kN�m)
        if a(i) < h_1
            M_real(i) = (P_con * (y_c_1 - a(i) / 2) - sum(ms_i)) * 0.0000980665;
        else
            M_real(i) = (P_con_ala * (y_c_1 - h_1 / 2) + P_con_alma * (y_c_1 - (a(i) - h_1) / 2) - sum(ms_i)) * 0.0000980665;
        end

        % Calculate the axial force (Pn) in kNs
        P_real(i) = (P_con - sum(f_s)) * 0.00980665;
        
        % Check if the axial force exceeds the maximum compressive load
        if P_n_max >= P_real(i)
            P_n(i) = P_real(i);
            M_n(i) = M_real(i);
        else
            M_n(i) = 0;
            P_n(i) = P_n_max;
        end
    end
    
    % Calculate the transverse forces and moments
    P_viga = 0.1 * fpc * Ag * 0.00980665;
    M_b = max(M_n);
    p_1 = find(M_n == M_b); % Position of the maximum moment
    P_b = P_n(p_1);
    
    P_tran_1 = P_viga / 0.7; % For walls with stirrups
    P_tran_2 = P_b / 0.7;
    
    P_tran = min(P_tran_1, P_tran_2);

    % Calculate phi_Mn and phi_Pn
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
    
    % Store data for the second half of the diagram
    flag = flag + 1;
    
    if flag == 0 || flag == 1
        M_n_neg = -M_n;
        P_n_neg = P_n;
        e_i_a_neg = e_i_a;
        
        phi_Mn_neg = -phi_Mn;
        phi_Pn_neg = phi_Pn;
        
        % Swap variables for the negative axis of the diagram
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
        M_n_pos = M_n;
        P_n_pos = P_n;
        e_i_a_pos = e_i_a;
        
        phi_Mn_pos = phi_Mn;
        phi_Pn_pos = phi_Pn;
    end
end

% Reverse the order of the graph data to close the interaction diagram
M_n_pos_graph = [];
P_n_pos_graph = [];
phi_M_n_pos_graph = [];
phi_P_n_pos_graph = [];
n_3 = length(M_n_pos);

for i_4 = 0:n_3-1
    M_n_pos_graph = [M_n_pos_graph; M_n_pos(n_3 - i_4)];
    P_n_pos_graph = [P_n_pos_graph; P_n_pos(n_3 - i_4)];
    phi_M_n_pos_graph = [phi_M_n_pos_graph; phi_Mn_pos(n_3 - i_4)];
    phi_P_n_pos_graph = [phi_P_n_pos_graph; phi_Pn_pos(n_3 - i_4)];
end

% Combine positive and negative values for the full interaction diagram
M_n_graph = [M_n_neg; M_n_pos_graph];
P_n_graph = [P_n_neg; P_n_pos_graph];

phi_M_n_graph = [phi_Mn_neg; phi_M_n_pos_graph];
phi_P_n_graph = [phi_Pn_neg; phi_P_n_pos_graph];

%% Plot and Save Results
% Nominal interaction diagram
figure(1)
plot(M_n_pos, P_n_pos, 'LineWidth', 2)
hold on
plot(M_n_neg, P_n_neg, 'LineWidth', 2)
xlabel('Moment (kN�m)')
ylabel('Axial Load (kN)')
title('Interaction Diagram - T-Beam')
grid on
legend('Positive Axis', 'Negative Axis')
close(gcf)

% Phi-factor modified diagram
figure(2)
plot(M_n_graph, P_n_graph, 'LineWidth', 2)
hold on
plot(phi_M_n_graph, phi_P_n_graph, '--', 'LineWidth', 2)
xlabel('Moment (kN�m)')
ylabel('Axial Load (kN)')
title('Interaction Diagram with phi Factor')
grid on
legend('Nominal Values', 'Modified Values with phi')
saveas(gcf, fullfile(outputDir, 'phi_modified_interaction.png'));
close(gcf)