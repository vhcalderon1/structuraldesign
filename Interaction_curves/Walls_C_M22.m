%% INTERACTION DIAGRAM GENERATOR FOR REINFORCED CONCRETE STRUCTURAL WALLS
% 
% Author: Victor Calderón (November 2019)
% Updated: Jefferson De la Cuba (February 2025)
% ========================================================================
%
% Generates nominal/factored interaction diagrams per ACI-318 requirements.
% Analyzes combined axial-flexural capacity considering material nonlinearity,
% reinforcement layout, and strain compatibility.
%
% INPUTS:
% - Material: fpc (kg/cm²), fy (kg/cm²), Es (kg/cm²), strains
% - Geometric: Wall segment dimensions (cm)
% - Reinforcement: Layer areas (cm²), spacing (cm)
% - Dataset: ../datasets/C_M22.txt (Applied loads: kN, kN-m)
%
% OUTPUTS:
% - Interaction diagrams (PNG images)
% - Console output of maximum axial capacity
% - MATLAB figure objects of P-M relationships
%
% ========================================================================

%% Initialize workspace
clear
close all
clc

%% Configure output directory
if exist('../outputs', 'dir') == 0
    mkdir('../outputs');
end
output_path = '../outputs/';

%% Input - Material Properties
fpc = 210;         % Concrete compressive strength (kg/cm^2)
beta_c = 0.85;     % Stress block factor for concrete compression
fy = 4200;         % Steel yield strength (kg/cm^2) 
Es = 2.0000e+06;   % Steel modulus of elasticity (kg/cm^2)
ey = 0.0020;       % Steel yield strain
ecu = 0.003;       % Ultimate concrete strain (per ACI-318)

%% Input - Geometric Properties (cm)
% Element 1 - Top rectangle (smaller dimension)
% Element 2 - Web
% Element 3 - Bottom rectangle (larger dimension)

b1 = 130;  % Width 1
b2 = 10;   % Web width
b3 = 172;  % Width 3
h1 = 10;   % Height 1
h2 = 285;  % Web height
h3 = 10;   % Height 3
Ag = b1*h1 + b2*h2 + h3*b3;  % Gross cross-sectional area

%% Input - Steel Distribution
n = 17;  % Number of reinforcement layers

% From bottom to top:
As = [repmat(1.13,7,1); 0.5; 1.13; 1.13; 0.5; 1.13; ones(3,1); 2.26; 9.76];  % Reinforcement areas (cm^2)

eid = [3; 12; 12; 12; 6; 6; 6; 1; 11; 12; 2; 10; 11; 25; 25; 8; 5];  % Spacing between bars (cm)
ei = [];  % Initialize bar positions (cm)

% Calculate positions along the section (cumulative spacing)
for i = 1:length(eid)-1
    ei(1) = eid(1);
    ei = [ei, ei(i) + eid(i+1)];
end
ei = ei';

n_c = 1;  % Increment for neutral axis depth

%% Compression Analysis (Negative Side)
% Calculate the centroid of the gross section (from top)
y_c2 = 0.5*(h1*b1^2 + h2*b2^2 + h3*b3^2) / (b1*h1 + b2*h2 + b3*h3);
y_c1 = b3 - y_c2;  % Distance from the bottom

% Define a range of neutral axis depths to evaluate
c = 0.00001:n_c:b3;  % (cm) – adjust as needed

e_i_a = [];  % Store strains for each neutral axis position

As_2 = [];
ei_2 = [];

% Reverse the order of reinforcement for the compression side (from top)
for i_3 = 0:length(As)-1
    As_2 = [As_2; As(n - i_3)];
    ei_2 = [ei_2; b3 - ei(n - i_3)];
end

flag = 0;

% Maximum compressive load with ties (kN)
P_n_max = 0.8*(0.85*fpc*(Ag - sum(As)) + fy*sum(As)) * 0.00981;
lim_1 = b3 - b1;
lim_2 = b3 - b2;

% Initialize arrays
a         = zeros(length(c),1);  % Depth of equivalent stress block
P_n       = zeros(length(c),1);  % Nominal axial force (kN)
P_real    = zeros(length(c),1);  % Actual axial force (kN)
M_n       = zeros(length(c),1);  % Nominal moment (kN-m)
M_real    = zeros(length(c),1);  % Actual moment (kN-m)
phi_Mn    = zeros(length(c),1);  % Factored moment (kN-m)
phi_Pn    = zeros(length(c),1);  % Factored axial force (kN)

for i = 1:length(c)
    e_i = zeros(length(As),1);
    f_s = zeros(length(As),1);
    ms_i = zeros(length(As),1);
    
    a(i) = c(i) * beta_c;  % Depth of the equivalent stress block
    
    for i_1 = 1:length(As)
        % Calculate strain in each reinforcement layer
        e_i(i_1) = ecu*(ei(i_1) - c(i)) / c(i);
        % (Negative: Compression, Positive: Tension)
        
        % Calculate force in each reinforcement layer (kg)
        if abs(e_i(i_1)) > ey
            if e_i(i_1) > 0
                f_s(i_1) = fy*As(i_1);
            else
                f_s(i_1) = -fy*As(i_1);
            end
        else
            f_s(i_1) = Es*e_i(i_1)*As(i_1);
        end
        
        % Calculate concrete compressive force (kg)
        if a(i) < lim_1
            P_con = 0.85*fpc*a(i)*h3;
        elseif lim_1 <= a(i) && a(i) < lim_2
            P_con_3 = 0.85*fpc*h3*a(i);
            P_con_1 = 0.85*fpc*(a(i) - (b3 - b1))*h1;
            P_con = P_con_3 + P_con_1;
        elseif a(i) >= lim_2
            P_con_1 = 0.85*fpc*h1*(a(i) - (b3 - b1));
            P_con_2 = 0.85*fpc*h2*(a(i) - (b3 - b2));
            P_con_3 = 0.85*fpc*h3*a(i);
            P_con = P_con_1 + P_con_2 + P_con_3;
        end
        
        % Moment contribution from each reinforcement layer (kg-cm)
        ms_i(i_1) = f_s(i_1) * (y_c1 - ei(i_1));
    end
    
    % Store the strains for this neutral axis depth
    e_i_a = [e_i_a; e_i'];
    
    % Calculate the nominal moment (M_n) in kN-m
    if a(i) < lim_1
        M_real(i) = (P_con*(y_c1 - a(i)/2) - sum(ms_i)) * 0.0000980665;
    elseif lim_1 <= a(i) && a(i) < lim_2
        M_real(i) = (P_con_1*(y_c1 - (a(i) + b3 - b1)/2) + P_con_3*(y_c1 - a(i)/2) - sum(ms_i)) * 0.0000980665;
    elseif a(i) >= lim_2
        M_real(i) = (P_con_1*(y_c1 - (a(i) + b3 - b1)/2) + P_con_2*(y_c1 - (a(i) + b3 - b2)/2) + P_con_3*(y_c1 - a(i)/2) - sum(ms_i)) * 0.0000980665;
    end
    
    % Total axial force (P_n) in kN (Positive: Compression, Negative: Tension)
    P_real(i) = (P_con - sum(f_s)) * 0.00981;
    
    if P_n_max >= P_real(i)
        P_n(i) = P_real(i);
        M_n(i) = M_real(i);
    else
        M_n(i) = 0;
        P_n(i) = P_n_max;
    end
end

% Beam (or wall) axial load (kN)
P_beam = 0.1 * fpc * Ag * 0.00981;
M_b = max(M_n);
p_1 = find(M_n == M_b);  % Index of maximum moment
P_b = P_n(p_1);

P_tran_1 = P_beam / 0.7;  % For walls with ties
P_tran_2 = P_b / 0.7;
P_tran = min(P_tran_1, P_tran_2);

% Calculate factored capacities (phi_Mn and phi_Pn)
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

% Store the negative side data (for the negative moment diagram)
M_n_neg = -M_n;
P_n_neg = P_n;
e_i_a_neg = e_i_a;

phi_Mn_neg = -phi_Mn;
phi_Pn_neg = phi_Pn;

%% Switch to the Positive Side of the Diagram
% Reassign reinforcement data for the positive (tension) side
As = As_2;
ei = ei_2;
lim_1 = b2;
lim_2 = b1;
y_c1 = y_c2;

% Reinitialize arrays for the positive side
a         = zeros(length(c),1);
P_n       = zeros(length(c),1);
P_real    = zeros(length(c),1);
M_n       = zeros(length(c),1);
M_real    = zeros(length(c),1);
phi_Mn    = zeros(length(c),1);
phi_Pn    = zeros(length(c),1);
e_i_a     = [];

for i = 1:length(c)
    e_i = zeros(length(As),1);
    f_s = zeros(length(As),1);
    ms_i = zeros(length(As),1);
    
    a(i) = c(i) * beta_c;
    
    for i_1 = 1:length(As)
        % Calculate strain in each reinforcement layer
        e_i(i_1) = ecu*(ei(i_1) - c(i)) / c(i);
        
        % Calculate force in each reinforcement layer (kg)
        if abs(e_i(i_1)) > ey
            if e_i(i_1) > 0
                f_s(i_1) = fy*As(i_1);
            else
                f_s(i_1) = -fy*As(i_1);
            end
        else
            f_s(i_1) = Es*e_i(i_1)*As(i_1);
        end
        
        % Calculate concrete compressive force (kg)
        if a(i) < lim_1
            P_con = 0.85*fpc*a(i)*(h1 + h2 + h3);
        elseif lim_1 <= a(i) && a(i) < lim_2
            P_con_2 = 0.85*fpc*h2*b2;
            P_con_1 = 0.85*fpc*a(i)*(h3 + h1);
            P_con = P_con_2 + P_con_1;
        elseif a(i) >= lim_2
            P_con_1 = 0.85*fpc*h1*b1;
            P_con_2 = 0.85*fpc*h2*b2;
            P_con_3 = 0.85*fpc*h3*a(i);
            P_con = P_con_1 + P_con_2 + P_con_3;
        end
        
        % Moment contribution from each reinforcement layer (kg-cm)
        ms_i(i_1) = f_s(i_1) * (y_c1 - ei(i_1));
    end
    
    % Store the strains
    e_i_a = [e_i_a; e_i'];
    
    % Calculate the nominal moment (M_n) in kN-m
    if a(i) < lim_1
        M_real(i) = (P_con*(y_c1 - a(i)/2) - sum(ms_i)) * 0.0000980665;
    elseif lim_1 <= a(i) && a(i) < lim_2
        M_real(i) = (P_con_1*(y_c1 - a(i)/2) + P_con_2*(y_c1 - b2/2) - sum(ms_i)) * 0.0000980665;
    elseif a(i) >= lim_2
        M_real(i) = (P_con_1*(y_c1 - b1/2) + P_con_2*(y_c1 - b2/2) + P_con_3*(y_c1 - a(i)/2) - sum(ms_i)) * 0.0000980665;
    end
    
    % Total axial force (kN)
    P_real(i) = (P_con - sum(f_s)) * 0.00981;
    
    if P_n_max >= P_real(i)
        P_n(i) = P_real(i);
        M_n(i) = M_real(i);
    else
        M_n(i) = 0;
        P_n(i) = P_n_max;
    end
end

% Beam (or wall) axial load (kN)
P_beam = 0.1 * fpc * Ag * 0.00981;
M_b = max(M_n);
p_1 = find(M_n == M_b);
P_b = P_n(p_1);

P_tran_1 = P_beam / 0.7;
P_tran_2 = P_b / 0.7;
P_tran = min(P_tran_1, P_tran_2);

% Calculate factored capacities (phi_Mn and phi_Pn) for the positive side
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

M_n_pos = M_n;
P_n_pos = P_n;
e_i_a_pos = e_i_a;

phi_Mn_pos = phi_Mn;
phi_Pn_pos = phi_Pn;

% Rearranging data for closing the interaction diagram (reverse the positive side)
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

% Combine negative and positive sides to form the complete interaction diagram
M_n_graf = [M_n_neg; M_n_pos_graf];
P_n_graf = [P_n_neg; P_n_pos_graf];

phi_M_n_graf = [phi_Mn_neg; phi_M_n_pos_graf];
phi_P_n_graf = [phi_Pn_neg; phi_P_n_pos_graf];

%% Plotting the Interaction Diagrams

% Figure 1: Positive/Negative Reinforcement Comparison
figure
plot(M_n_pos, P_n_pos, 'b-', 'LineWidth', 2)
hold on
plot(M_n_neg, P_n_neg, 'r-', 'LineWidth', 2)
grid on
xlabel('Moment (kN-m)')
ylabel('Axial Load (kN)')
title('Positive vs. Negative Reinforcement Interaction')
legend('Positive', 'Negative', 'Location', 'best')
close(gcf)

% Figure 2: Nominal vs Factored Capacity
figure
plot(M_n_graf, P_n_graf, 'Color', [0 0.8 1], 'LineWidth', 2)
hold on
plot(phi_M_n_graf, phi_P_n_graf, '--', 'Color', [0 0.8 1], 'LineWidth', 2)
grid on
xlabel('Moment (kN-m)')
ylabel('Axial Load (kN)')
title('Nominal vs Factored Interaction Diagram')
legend('Nominal', 'Factored', 'Location', 'best')
saveas(gcf, fullfile(output_path, 'Figure2_NominalFactored.png'))
close(gcf)

%% Applied Loads on the Wall

% Load applied loads data from the external file (relative path to datasets folder)
load('../datasets/C_M22.txt')
P_u_M22 = C_M22(:,1);  % Applied axial loads (kN)
M_u_M22 = C_M22(:,2);  % Applied moments (kN-m)
sz = 15;             % Marker size for scatter plot

% Define plot limits
x_lim    = 0.5 * 1e3;  % kN-m
y_lim_sup = 2 * 1e3;   % kN
y_lim_inf = 0.5 * 1e3; % kN

% Figure 3: Applied Loads Visualization
figure
plot(M_n_graf, P_n_graf, 'Color', [0 0.8 1], 'LineWidth', 2)
hold on
plot(phi_M_n_graf, phi_P_n_graf, '--', 'Color', [0 0.8 1], 'LineWidth', 2)
scatter(M_u_M22, P_u_M22, 15, [0.7 0.7 0.7], 'filled')
plot([-x_lim, x_lim], [0, 0], 'k-', 'LineWidth', 0.9)
plot([0, 0], [-y_lim_inf, y_lim_sup], 'k-', 'LineWidth', 0.9)
grid on
xlabel('Moment (kN-m)')
ylabel('Axial Load (kN)')
title('Interaction Diagram with Applied Loads')
legend('Nominal', 'Factored', 'Applied Loads', 'Location', 'best')
saveas(gcf, fullfile(output_path, 'Figure3_AppliedLoads.png'))
close(gcf)

%% Calculation of Maximum Axial Load
Pu_Nch = Ag * fpc * 0.3 * 0.00981;  % Maximum axial load (kN)