%% Pressure Drop in DX-SAHP (Based on R134a)
%clear; clc;
%Installing the library for coolprop or something idk
%[v,e] = pyversion; system([e, ' -m pip install --user -U CoolProp']);

% Lengths of Sections
num_sec = 4; % Number of sections
L = [1 1 1 1]; % Lengths between components [m]
% L1 = Length between 1 (Compressor) and 2 (Condenser Inlet)
% L2 = Length between 2 (Condenser Exit) and 3 (Expansion Valve)
% L3 = Length between 3 (Expansion valve) and 4 (Evaporator Inlet)
% L4 = Length between 4 (Evaporator Exit) and 1 (Compressor)

% Pipe Parameters
D1 = 1; %inner pipe diameter [in]
D_pipe = [D1 D1 D1 D1]*0.0254; % Inside Pipe Diameters [in]*0.0254 -> [m]
k = 1.524E-6; % Roughness coefficient of drawn copper piping [m]

% Refrigerant Parameters
T = [333 283 271 295]; % Temperature [K]
Q = [1 0 1 1]; % Quality **quality of mixture at 3 is unknown and coolprop didn't accept two-phase**
rho = zeros(1,num_sec); % Density Matrix [kg/m^3]
mu = zeros(1,num_sec); % Viscosity Matrix [Pas]
v = zeros(1,num_sec); % Kinematic Viscosity Matrix [m^2/s]
% Pressures 8E5 Pa and 3.76E6 Pa

% Looping for each section
for i = 1:num_sec
    rho(i) = py.CoolProp.CoolProp.PropsSI('D', 'T', T(i), 'Q', Q(i), 'R134a'); % Density [kg/m^3]
    mu(i) = py.CoolProp.CoolProp.PropsSI('V', 'T', T(i), 'Q', Q(i), 'R134a'); % Viscosity [Pas]
    v(i) = mu(i)/rho(i); % Kinematic Viscosity [m^2/s]
end

%Flow Parameters
m_dot = ones(1, num_sec)*0.04; % Mass flow rate [kg/s]
V_dot = m_dot./rho; % Volumetric flow rate [m3/s]
w = V_dot./(pi*(D_pipe/2).^2); % Velocity of fluid [m/s]
Re = w.*D_pipe./v; % Reynold's Number
% Re_Test = rho.*w.*D_pipe./mu; % Testing reynold's number
Laminar = zeros(1, num_sec);

% Determining Flow Type
for i = 1:num_sec
    if Re(i) < 2320
        Laminar(i) = 1;
    else
        Laminar(i) = 0;
    end
end

% Pressure Drops
f_coeff = zeros(1, num_sec); % Pipe friction coefficient
deltaP_fric = zeros(1, num_sec); % Pressure drop due to friction [Pa]
deltaP_comps_T = zeros(1, num_sec); % Pressure drop from components [Pa]
deltaP_comps_L = zeros(1, num_sec);
deltaP_comps = zeros(1, num_sec);
deltaP_height = zeros(1, num_sec); % Pressure drop from height [Pa]

% Calculating Pressure drop due to friction
for i = 1:num_sec
    if Laminar(i) == 1 % Laminar Flow
        f_coeff(i) = 64/Re(i); % Pipe friction coefficient
        deltaP_fric(i) = (f_coeff(i) * L(i) * rho(i) * w(i)^2) / (D_pipe(i)*2); % Pressure Drop [Pa}
    elseif Laminar == 0 % Turbulent Flow
        f_coeff(i) = fzero( @(f) 1/sqrt(f) + 2*log10(((k/D_pipe(i))/3.7) + (2.51/(Re(i)*sqrt(f)))), [1E-18, 1]);
        deltaP_fric(i) = (f_coeff(i) * L(i) * rho(i) * w(i)^2) / (D_pipe(i)*2); % [Pa]
    end
end

% Pressure drop due to bends and components
num_90Lbends = [1 1 1 1]; % Number of 90° Long bends for each section [1 1 1 1]
num_90Tbends = [1 1 1 1]; % Number of 90° threaded bends for each section [0 0 0 0]
K_90_L = 0.2; % Loss coefficient for a long radius flanged 90° elbow
K_90_T = 1.5; % Loss coefficient for a threaded 90° elbow

% Pressure drop for all components per section [Pa]
for i = 1:num_sec
    deltaP_comps_T(i) = (num_90Tbends(i) * K_90_T * rho(i) * w(i)^2) / 2;
    deltaP_comps_L(i) = (num_90Lbends(i) * K_90_L * rho(i) * w(i)^2) / 2;
    deltaP_comps(i) = deltaP_comps_T(i) + deltaP_comps_L(i);
    
end

% Pressure drop/gain due to height
height_diff = [0 0 0 0]; % Relative difference in height between components [m]
g = 9.81; % Gravity [m/s^2]

% Pressure drop per section for change in height [Pa]
for i = 1:num_sec
    deltaP_height(i) = rho(i) * g * height_diff(i);
end

% Total Pressure drop for each section [Pa]
deltaP_Total = deltaP_comps + deltaP_fric + deltaP_height;

%% Pressure Drop From Friction Plots
sections = [1 2 3 4];

% bar(sections, deltaP_Total_quarter_in, 1); hold on
% bar(sections, deltaP_Total_half_in, 1);
% bar(sections, deltaP_Total_3_4_in, 1);
% bar(sections, deltaP_Total_one_in, 1);
P_drop_combined = [deltaP_fric_25
                   deltaP_fric_5_16
                   deltaP_fric_3_8
                   deltaP_fric_50
                   deltaP_fric_75
                   deltaP_fric_1];
P_drop_combined = P_drop_combined';
bar(sections,P_drop_combined, 'grouped');
legend("1/4 inch", "5/16 inch", "3/8 inch", "1/2 inch", "3/4 inch", "1 inch");
xlabel("Section");
ylabel("Pressure Drop, [Pa/m]");
%title("Pressure Drop Due to Friction Per Meter of Piping For 1/2, 3/4, and 1 inch Pipe Diameters");

%% Pressure Drop From bends Plots

sections = [1 2 3 4];
P_drop_combined = [deltaP_comps_1_90L
                   deltaP_comps_1_90T
                  ];
              
P_drop_combined = P_drop_combined';
bar(sections,P_drop_combined, 'grouped');
legend("Long 90° Elbow", "Threaded 90° Elbow");
xlabel("Section");
ylabel("Pressure Drop, [Pa]");
%title("Pressure Drop Due to a Long 90° Elbow & Threaded 90° Elbow for a 1 inch Pipe");





