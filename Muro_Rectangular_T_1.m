%% T-WALL INTERACTION DIAGRAM GENERATOR
% 
% Author: Victor Calderón (July 2019)
% Updated: Jefferson De la Cuba (February 2025)
% --------------------------------------------------------------------------
% Generates nominal/design interaction diagrams for reinforced concrete T-walls per ACI 318
% Includes validation against applied load cases and exports results to specified directory
%
% Inputs:
%   - Material: fpc (kg/cm²), fy (kg/cm²), Es (kg/cm²), strains (ey, ecu)
%   - Geometry: b1/b2 (cm), h1/h2 (cm), steel layer coordinates (ei)
%   - Reinforcement: As (cm²) for each layer, bar spacing (eid)
%   - Dataset: Rec_P1_M33.txt (contains applied P [kN] and M [kN-m])
%
% Outputs:
%   - MATLAB figures: Interaction diagrams in .fig format
%   - JPEG image: Final annotated diagram in outputs directory
%   - Console output: Maximum axial load capacity (kN)

clear
close all
clc

%% Directory Configuration
% Create outputs directory if non-existent (one level up from current file)

outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'outputs'); % Robust path traversal
if ~exist(outputDir, 'dir')
    mkdir(outputDir); 
end

% Dataset path configuration (outside code directory)
dataFile = fullfile('..', 'datasets', 'Rec_P1_M33.txt');

%% Interaction Diagram - T Wall

% Input - Material Properties
fpc    = 210;           % Concrete compressive strength (kg/cm²)
beta_c = 0.85;          % Compression block factor for concrete
fy     = 4200;          % Steel yield strength (kg/cm²)
Es     = 2.0000e+06;    % Modulus of elasticity of steel (kg/cm²)
ey     = 0.0020;        % Steel yield strain
ecu    = 0.003;         % Ultimate strain of concrete (per ACI-318)

% Conversion factors
kgf_to_kN = 0.00980665;    % 1 kgf = 0.00980665 kN
kgfcm_to_kNm = 0.0000980665; % 1 kgf·cm = 0.0000980665 kN·m

% Input - Geometric Properties (cm)
b_1 = 0;                % Flange width
b_2 = 15;               % Web width
h_1 = 0;                % Flange height
h_2 = 675;              % Total height (web height = h_2 - h_1)
Ag  = b_1*h_1 + b_2*(h_2-h_1);   % Gross area (cm²)

%% Input - Distribution of Steel Reinforcement
n = 41;  % Number of steel layers

% From top to bottom:
As = [repmat(2.26,5,1); ones(31,1); repmat(2.26,5,1)];  % Steel area (cm²)
eid = [3; repmat(7,5,1); repmat(20,30,1); repmat(7,5,1)];  % Spacing between bars (cm)

% Distance from top to each steel layer (cm)
ei = zeros(length(eid),1);
ei(1) = eid(1);
for i = 2:length(eid)
    ei(i) = ei(i-1) + eid(i);
end

n_c = 1;  % Step increment for neutral axis depth

%% Compression Analysis
% Centroid of gross section (from top)
y_c_1 = (0.5*((b_1*h_1^2) + (h_2^2 - h_1^2)*b_2)) / (b_1*h_1 + b_2*(h_2-h_1));
y_c_2 = h_2 - y_c_1;

% Neutral axis depth range (cm)
c = 0.00001 : n_c : h_2 + 0.5*h_2;  

% Reverse reinforcement for negative bending
As_2 = As(end:-1:1);
ei_2 = h_2 - ei(end:-1:1);

flag = 0;

% Maximum compressive load (kN)
P_n_max = 0.8*(0.85*fpc*(Ag - sum(As)) + fy*sum(As)) * kgf_to_kN;

while flag <= 1
    a       = zeros(length(c),1);
    P_n     = zeros(length(c),1);
    P_real  = zeros(length(c),1);
    M_n     = zeros(length(c),1);
    M_real  = zeros(length(c),1);
    phi_Mn  = zeros(length(c),1);
    phi_Pn  = zeros(length(c),1);

    for i = 1:length(c)
        e_i  = zeros(length(As),1);
        f_s  = zeros(length(As),1);
        ms_i = zeros(length(As),1);

        a(i) = c(i) * beta_c;  % Depth of stress block

        % Steel layer forces
        for i_1 = 1:length(As)
            e_i(i_1) = ecu*(ei(i_1) - c(i)) / c(i); % Strain
            
            if e_i(i_1) > ey
                f_s(i_1) = fy * As(i_1); % Tension yield
            elseif e_i(i_1) < -ey
                f_s(i_1) = -fy * As(i_1); % Compression yield
            else
                f_s(i_1) = Es * e_i(i_1) * As(i_1); % Elastic
            end
            
            ms_i(i_1) = f_s(i_1) * (y_c_1 - ei(i_1)); % Moment contribution
        end

        % Concrete force
        if a(i) < h_1
            P_con = 0.85 * fpc * a(i) * b_1;
        else
            P_con_ala  = 0.85 * fpc * h_1 * b_1;
            P_con_alma = 0.85 * fpc * (a(i)-h_1) * b_2;
            P_con = P_con_ala + P_con_alma;
        end

        % Convert forces and moments
        P_real(i) = (P_con - sum(f_s)) * kgf_to_kN; % kN
        if a(i) < h_1
            M_real(i) = (P_con * (y_c_1 - a(i)/2) - sum(ms_i)) * kgfcm_to_kNm;
        else
            M_real(i) = (P_con_ala*(y_c_1 - h_1/2) + P_con_alma*(y_c_1 - (a(i)-h_1)/2) - sum(ms_i)) * kgfcm_to_kNm;
        end
        
        % Enforce Pn_max
        if P_n_max >= P_real(i)
            P_n(i) = P_real(i);
            M_n(i) = M_real(i);
        else
            M_n(i) = 0;
            P_n(i) = P_n_max;
        end
    end
    
    % Transition load calculations
    P_viga = 0.1 * fpc * Ag * kgf_to_kN;  % kN
    M_b    = max(M_n);
    p_1    = find(M_n == M_b, 1);
    P_b    = P_n(p_1);
    P_tran = min(P_viga/0.7, P_b/0.7);
    
    % Strength reduction factors
    for i_5 = 1:length(c)
        phi = 0.9 - 0.2 * P_n(i_5)/P_tran;
        phi = max(0.7, min(0.9, phi));
        phi_Mn(i_5) = M_n(i_5) * phi;
        phi_Pn(i_5) = P_n(i_5) * phi;
    end
    
    % Swap parameters for negative side
    if flag == 0
        M_n_neg      = -M_n;
        P_n_neg      = P_n;
        phi_Mn_neg   = -phi_Mn;
        phi_Pn_neg   = phi_Pn;
        As    = As_2;
        ei    = ei_2;
        [h_1, h_2] = deal(h_2, h_1);
        [b_1, b_2] = deal(b_2, b_1);
        y_c_1 = y_c_2;
    else
        M_n_pos     = M_n;
        P_n_pos     = P_n;
        phi_Mn_pos  = phi_Mn;
        phi_Pn_pos  = phi_Pn;
    end
    flag = flag + 1;
end

% Combine interaction diagram results
M_n_graf     = [M_n_neg; flipud(M_n_pos)];
P_n_graf     = [P_n_neg; flipud(P_n_pos)];
phi_M_n_graf = [phi_Mn_neg; flipud(phi_Mn_pos)];
phi_P_n_graf = [phi_Pn_neg; flipud(phi_Pn_pos)];

%% Plotting
figure;
plot(M_n_graf, P_n_graf, 'b-', phi_M_n_graf, phi_P_n_graf, 'r--', 'LineWidth', 2);
grid on;
xlabel('Moment (kN·m)');
ylabel('Axial Load (kN)');
legend('Nominal', 'Design', 'Location', 'Best');
title('Interaction Diagram');

% Save as FIG
saveas(gcf, fullfile(outputDir, 'InteractionDiagram.png')); 

% Export PNG
print(fullfile(outputDir, 'C_I_Conven_TRec_M33.png'), '-dpng', '-r300'); 

close(gcf); % Always close after saving

%% Maximum Axial Load
Pu_Nch = Ag * fpc * 0.3 * kgf_to_kN;  % kN