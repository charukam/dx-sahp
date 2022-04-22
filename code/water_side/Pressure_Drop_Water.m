% This code calculates the pressure drop in a circular pipe
% clear;
% Installing the library for coolprop
% [v,e] = pyversion; system([e, ' -m pip install --user -U CoolProp']);
clc;

% Pipe Parameters
L_ft = 15; % [ft]
L = L_ft/3.281; % Length of piping [ft]/3.281 = [m]
D_pipe_in = 0.75; % [in]
D_pipe = D_pipe_in*0.0254; % Inside Pipe Diameters [in]*0.0254 -> [m]
k_copper = 1.524E-6; % Roughness coefficient of drawn copper piping [m]
% roughness coefficient values: https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html

% Initializing Fluid Parameters
T = 283.15; % Temperature [K]
Q = 0; % Quality (0 for saturated liquid)

% Obtaining density,vicosity, and kinematic viscosity values from coolprop 
rho = py.CoolProp.CoolProp.PropsSI('D', 'T', T, 'Q', Q, 'water'); % Density [kg/m^3]
mu = py.CoolProp.CoolProp.PropsSI('V', 'T', T, 'Q', Q, 'water'); % Viscosity [Pas]
v = mu/rho; % Kinematic Viscosity [m^2/s]

%Flow Parameters
m_dot = 0.5*0.063; % Mass flow rate [USgpm]*0.063 = [kg/s]
V_dot = m_dot/rho; % Volumetric flow rate [m3/s]
w = V_dot/(pi*(D_pipe/2)^2); % Velocity of fluid [m/s]
Re = w*D_pipe/v; % Reynold's Number
Re_Test = rho.*w.*D_pipe./mu; % Verifying reynold's number
Laminar = 0;

% Determining Flow Type
if Re < 2320
    Laminar = 1;
else
    Laminar = 0;
end

% Initializing Pressure Drop Variables
f_coeff = 0; % Pipe friction coefficient
deltaP_fric = 0; % Pressure drop due to friction [Pa]
deltaP_comps_T = 0; % Pressure drop from threaded bends [Pa]
deltaP_comps_L = 0; % Pressure drop from long radius bends [Pa]
deltaP_comps = 0; % Total Pressure drop from bends
deltaP_height = 0; % Pressure drop from height [Pa]

% Calculating Pressure drop due to friction
if Laminar == 1 % Laminar Flow
    f_coeff = 64/Re; % Circular pipe friction coefficient
    deltaP_fric = (f_coeff * L * rho * w^2) / (D_pipe*2); % Pressure Drop - Darcy-Weisbach [Pa]
elseif Laminar == 0 % Turbulent Flow - Colebrook White friction coefficient
    f_coeff = fzero( @(f) 1/sqrt(f) + 2*log10(((k_copper/D_pipe)/3.7) + (2.51/(Re*sqrt(f)))), [1E-18, 1]);
    deltaP_fric = (f_coeff * L * rho * w^2) / (D_pipe*2); % Pressure Drop - Darcy-Weisbach [Pa]
end

% Pressure drop due to bends and components
num_90Lbends = 0; % Number of 90° Long bends
num_90Tbends = 12; % Number of 90° threaded bends
num_bends_T = num_90Lbends + num_90Tbends; % Total number of bends
K_90_L = 0.2; % Loss coefficient for a long radius flanged 90° elbow
K_90_T = 1.5; % Loss coefficient for a threaded 90° elbow

% Pressure drop for all bends [Pa]
deltaP_comps_T = (num_90Tbends * K_90_T * rho * w^2) / 2;
deltaP_comps_L = (num_90Lbends * K_90_L * rho * w^2) / 2;
deltaP_comps = deltaP_comps_T + deltaP_comps_L;
    
% Pressure drop/gain due to height
height_diff_ft = 1.2; % [ft]
height_diff = height_diff_ft/3.281; % Relative difference in height between components [ft]/3.281 = [m]
g = 9.81; % Gravity [m/s^2]

% Pressure drop per section for change in height [Pa]
deltaP_height = rho * g * height_diff;

% Pressure drop from condenser coil
deltaP_Condenser_psi = 0.28; % Estimate from condenser coil data sheet
deltaP_Condenser_Pa = deltaP_Condenser_psi*6894.75729; % [psi]*6894.75729 = [Pa]

% Total Pressure drop [Pa]
deltaP_Total_Pa = deltaP_comps + deltaP_fric + deltaP_height + deltaP_Condenser_Pa;
deltaP_Total_psi = deltaP_Total_Pa/6894.75729;
deltaP_Total_ft = deltaP_Total_psi*2.31;

% Output
fprintf('Pressure drop for a %.2fin x %.fft pipe from:\n %.fft friction: %.2f Pa \n %.2fft height: %.2f Pa \n %.fx90° elbows: %.2f Pa \n Condenser: %.2f Pa \n'...
        ,D_pipe_in, L_ft, L_ft, deltaP_fric, height_diff_ft, deltaP_height, num_bends_T, deltaP_comps, deltaP_Condenser_Pa);
fprintf('\nThe total Pressure drop is %.2fPa = %.2fpsi = %.2fft \n',deltaP_Total_Pa,deltaP_Total_psi,deltaP_Total_ft);

