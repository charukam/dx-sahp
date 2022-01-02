clear
clc

%% Design Conditions
% Enviorment Data
import = readtable('winter.xlsx');
import = table2array(import);
T_a = import(:, 2);  % ambient temperature (K) at each irradiance point
I = import(:, 3); %instantenous irradiance (W/m^2)
T_1=T_a-2; %K
Twater_inlet=283.15; %K
Twater_outlet=328.15; %K
Water_consumtpion= 225; %Litres
Water_Density=1000; %g/L
Cp_Water= 4.18; %J/g.c
AvgOperational_hours=4.68; %hours
design_heatload=(Water_consumtpion.*Water_Density.*Cp_Water.*(Twater_outlet-Twater_inlet))./1000; %KJ
design_power=design_heatload./(4.68.*3600); %KW
T_conddesign=331.15; %k
ref='R410a';  


%% Thermodynamic Properties

%State 1 
for i= 1:length(T_1) 
H_1(i)=py.CoolProp.CoolProp.PropsSI('H','T',T_1(i),'Q',1, ref)/1000; %kJ/kg
S_1(i)=py.CoolProp.CoolProp.PropsSI('S','T',T_1(i),'Q',1, ref); %J/kgk
P_1(i)=py.CoolProp.CoolProp.PropsSI('P','T',T_1(i),'Q',1, ref); %Pa 
end 

H_1 = H_1';

%State 2
S_2=S_1; % Isentropic compression
P_2=py.CoolProp.CoolProp.PropsSI('P','T',T_conddesign,'Q',1,ref); %Pa 

for i=1:length(S_2) 

    H_2(i)=py.CoolProp.CoolProp.PropsSI('H','S',S_2(i),'P',P_2,ref)/1000; %KJ/kg
    T_2(i)=py.CoolProp.CoolProp.PropsSI('T','S',S_2(i),'P',P_2,ref); %K
end

H_2 = H_2';

%State 3
P_3=P_2; %Pa
T_3=py.CoolProp.CoolProp.PropsSI('T','P',P_3,'Q',0,ref); %K
H_3=py.CoolProp.CoolProp.PropsSI('H','P',P_3,'Q',0,ref)/1000; %J/kg


%State 4 
P_4=P_1;
H_4=H_3;


%% Thermodynamic Cycle Analysis

W_Compressor= H_2-H_1; %KJ/kg
Q_H=H_2-H_3; %KJ/Kg
Q_L=H_1-H_4; %KJ/Kg
COP=(Q_H)./(W_Compressor);

%% Solar Collector 
v_w = 14.2*1000/3600; %wind speed (m/s)
sigma = 5.67E-8; % STEFAN Boltzman Constant (m^2*kg/s^2*K)

%Collector Data
L_p=2; % length of plate (m)
W_p=1; %Width of plate (m)
A=L_p*W_p; %m^2
epsilon_c = 0.09; %emissivity of abosrber plate
T_fi = T_1; %evaporation temperatures (K) 
T_pm_check = T_fi+5; %first guess plate mean temperature (K)
alpha = 0.9; %absorbance of absorber
tau = 0.9; %transmittance of glazing
epsilon_g = 0.88; %emissivity of glazing material
T_pm  = 0;
tol = 0.1;
iterate = 0;

for i=1:length(T_fi) 
    Cp_ref_tfi(i) = py.CoolProp.CoolProp.PropsSI('C','T',T_fi(i),'Q',1, ref); %constant pressure specific heat of refridgerant at T_fi (J/kgK)
end
    Cp_ref_tfi = Cp_ref_tfi';

%Iteration to find plate mean temp
while (abs(T_pm_check - T_pm) > tol ) 
    
    T_pm = T_pm_check;    
    %Bottom Heat loss Coefficient
    delta_b = 0.05; % thickness of insulation (m)
    k_b = 0.035; %thermal conductivity of insulator (W/mK)
    U_b = k_b/delta_b; %bottom heat loss coefficient (W/m^2K)
    
    %Edge Heat Loss Coefficient
    k_e = 0.035; %thermal conductivity of insulator (W/mK)
    l_e = 0.05; %edge insulation thickness (m)
    delta_p= 0.004; %thickness of absorber plate (m)
    
    A_p = (2*L_p+2*W_p)*delta_p; %perimeter area (m^2)(2L+2w*Collector plate thickness)
    U_e = (k_e./l_e)*A_p./A; %edge heat loss coefficient (W/m^2K)
    
    %Top Heat Loss Coefficient
    h_w = 2.8 + 3*v_w; %convection due to wind (W/m^2*K)
    M = 1; %for single glazed collector
    beta = 45; %collector tilt angle (degrees)
    
    C = 520*(1-0.000051*beta^2);
    e = 0.43*(1 - 100./T_pm);
    f = (1+0.089*h_w - 0.116*h_w*epsilon_c)*(1+0.07866*M);
    
    U_tc = (M./((C./T_pm).*(((T_pm - T_a)/(M+f)).^e)) + 1/h_w).^-1; % heat loss due to convection( W/m^2.c)
    
    U_tr = (sigma*(T_pm.^2 + T_a.^2).*(T_pm + T_a))/((epsilon_c + 0.00591*M*h_w)^-1 + ((2*M + f - 1 +0.133*epsilon_c)/epsilon_g) - M); %heat loss due to radiation(( W/m^2.c)
    
    U_t = U_tc + U_tr; %total (W/m^2K)
    
    %Total Heat Loss
    U_L = U_e + U_b + U_t; %(W/m^2K)
    
    % Useful Heat Gain
    OD_tubes = 19.05/1000; %m
    tube_thickness = 0.8128/1000; % m for an M type copper tube 
    ID_tubes = OD_tubes - tube_thickness;
    Tube_Spacing = 100/1000; %m
    W = OD_tubes + Tube_Spacing; %m center to center distance between tubes 
    K_p = 205; %W/mK Thermal conductivity of aluminum absorber plate
    m = sqrt(U_L./(K_p *delta_p)); %for convenience
    h_fi = 300; %W/m^2.K %heat transfer coefficient between tube wall and refrigerant (unsure, should be between 300-600) 
    inv_cb = 0; % bond thermal conductivity assumed infinite
    F = tanh(m.*(W - OD_tubes)/2)./(m.*(W - OD_tubes)/2); %fin efficiency (unitless)
    F_prime_num = 1./U_L;
    F_prime_denom1 = 1./(U_L.*(OD_tubes + (W - OD_tubes).*F));
    F_prime_denom2 = inv_cb;
    F_prime_denom3 = 1/(pi*ID_tubes*h_fi);
    F_prime_denom = W.*(F_prime_denom1 + F_prime_denom2 + F_prime_denom3);
    F_prime = F_prime_num./F_prime_denom; %collector efficiency factor (unitless)
    Q_u = F_prime.*A.*(tau*alpha*I - U_L.*(T_fi - T_a)); %Hourly Useful Energy Gain (W)
    m_dot = (Q_u/1000)./Q_L; %mass flow rate (kg/s)
    
    %Check temperature
    FR = (m_dot.*Cp_ref_tfi)./(A.*U_L).*(1-exp(-A*U_L.*F_prime./(m_dot.*Cp_ref_tfi)));
    T_pm_check = T_fi + ((Q_u./A)./(FR.*U_L)).*(1 - FR);

    if ~isreal(T_pm_check)
        break
  

    end
iterate = iterate +1;
end

%% Efficiency
%Solar Collector Efficiency 

Q_u = F_prime.*A.*(tau*alpha*I - U_L.*(T_fi - T_a));
Q_h = A.*(tau*alpha*I - U_L.*(T_pm - T_a));
Collector_Efficiency = (Q_u)./(A*(I));
 
%Overall System Efficiency 
Overall_COP=Collector_Efficiency.*COP;

%Plot of System
figure
plot(T_1,Overall_COP);
figure 
plot(T_1,COP);
hold on
plot(T_1,Collector_Efficiency);


%% System properites 
RefMass_Flowrate= design_power./Q_H; %Kg/s 
Compressor_power=W_Compressor.*RefMass_Flowrate;
net_collectedheat_evap=(Q_L.*RefMass_Flowrate)./Collector_Efficiency;
