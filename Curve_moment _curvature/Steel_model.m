% ===============================================================================
% Steel Stress-Strain Curve Generator: Bilinear Elastic-Perfectly Plastic Model
% ==============================================================================
% Author: Victor Calderón (November 2020)
% Updated: Jefferson De la Cuba (February 2025)
% -------------------------------------------------------------------------
% Generates a bilinear elastic-perfectly plastic stress-strain curve for structural steel 
% using SI units. Implements Hookean elasticity with perfect plasticity post-yield. 
% Auto-configures output directory
%
% Input Parameters:
%   - yield_stress_MPa    : Yield stress [MPa] (scalar, typical range: 250-550 MPa)
%   - young_modulus_GPa   : Elastic modulus [GPa] (scalar, 200 GPa for steel)
%   - ultimate_strain     : Failure strain (scalar, typically 0.05-0.15)
%   - material_params.txt : Text file containing material parameters in the following order:
%                           1. Yield stress (MPa)        (e.g., 420)
%                           2. Young's modulus (GPa)     (e.g., 200)
%                           3. Ultimate strain           (e.g., 0.055)
%                           Parameters must be numeric values on separate lines.
%   - strain_data.txt     : (Optional) Custom strain input file for curve resolution control.
%
% Outputs:
%   - PNG figure          : Saved in '../outputs/' directory
%   - MATLAB figure       : Interactive plot window with dimensionless strain vs MPa stress
% -------------------------------------------------------------------------

%% Initialize Environment
clear
close all
clc

%% ====================== MATERIAL PARAMETERS =============================
yield_stress_MPa    = 420;       % ASTM A36 steel yield stress
young_modulus_GPa   = 200;       % Typical structural steel modulus
ultimate_strain     = 0.055;     % Ductility parameter

%% ====================== UNIT CONVERSIONS ================================
young_modulus_MPa = young_modulus_GPa * 1e3;  % GPa → MPa

%% ====================== STRAIN CALCULATIONS =============================
yield_strain = yield_stress_MPa / young_modulus_MPa;

elastic_strain = linspace(0, yield_strain, 100)';
plastic_strain = linspace(yield_strain, ultimate_strain, 100)';

%% ====================== STRESS CALCULATIONS =============================
elastic_stress = young_modulus_MPa * elastic_strain;
plastic_stress = repmat(yield_stress_MPa, 100, 1);

stress_curve = [elastic_stress; plastic_stress];
strain_curve = [elastic_strain; plastic_strain];

%% ====================== PLOT CONFIGURATION ==============================
figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 5 12 10])
plot(strain_curve, stress_curve, 'k-', 'LineWidth', 2)
grid on

% Axis labels with SI units
xlabel('Engineering Strain, $\varepsilon$',...
    'FontSize', 13,...
    'FontName', 'Times New Roman',...
    'Interpreter', 'latex')
ylabel('Engineering Stress, $\sigma$ (MPa)',...
    'FontSize', 13,...
    'FontName', 'Times New Roman',...
    'Interpreter', 'latex')

axis([-0.0025 0.06 -50 1.1*yield_stress_MPa])

% Legend correction (single entry)
legend('Elastic-Perfectly Plastic Model',...
    'Location', 'southeast',...
    'FontName', 'Times New Roman',...
    'Box', 'off')

%% ====================== FIGURE EXPORT ===================================
outputDir = fullfile('..', 'outputs');
if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

print(fullfile(outputDir, 'SteelStressStrain_Model.png'),...
    '-dpng',...
    '-r600')
close(gcf)