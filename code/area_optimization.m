clear
clc

n = 2; %number of transient irradiance data points

%Enviorment Data
T_a =  zeros(1, n);  % ambient temperature (K) at each irradiance point
I = zeros(1, n); %instantenous irradiance (W/m^2)
v_w = 10*1000/3600; %wind speed (m/s)
sigma = 1.38064852E-23; %Boltzman Constant (m^2*kg/s^2*K)

%Collector Data
A = 1:0.1:3; %m^2
epsilon_c = 0.1; %emissivity of abosrber plate
T_c = zeros(1, n); % mean temperature (K) of collector at each irradiance point
alpha = 0.08; %absorbance of glazing
tau = 1; %transmittance of glazing
epsilon_g = 0.92; %emissivity of glazing material

%Bottom Heat loss Coefficient
delta_b = 0.2; % thickness of insulation (m)
k_b = 0.1; %thermal conductivity of insulator ( )
h_b = 13; %convection coeff btw bottom and enviorment (W/m^2*K) between 12-25
U_b = delta_b/k_b +1/h_b; %bottom heat loss coefficient

%Edge Heat Loss Coefficient
h_e = 0.5; %convection coeff btw edge and enviorment (W/m^2*K)
A_p = 0.1; %perimeter area (m^2)(should be function of area, length, and width somehow)
U_e = h_e*A_p./A; %edge heat loss coefficient

%Top Heat Loss Coefficient
h_w = 5.7 + 3.8*v_w; %convection due to wind (W/m^2*K)
M = 1; %for single glazed collector
beta = 45; %collector tilt angle (degrees)

C = 520*(1-0.000051*beta^2);
e = 0.43*(1 - 100./T_c);
f = (1+0.089*h_w - 0.116*h_w*epsilon_c)*(1+0.07866*M);

U_tc = (M./((C./T_c).*((T_c - T_a)/(M+f))).^(exp(1)) + 1/h_w).^-1; % heat loss due to convection

U_tr = (sigma*(T_c.^2 + T_a.^2).*(T_c + T_a))/((epsilon_c + 0.059*M*h_w)^-1 + ((2*M - f - 1 +0.133*epsilon_c)/epsilon_g) - M); %heat loss due to radiation

U_t = sum(U_tc) + sum(U_tr); %total

%Total Heat Loss
U_L = U_e + U_b + U_t;
Q_o = U_L.*A*mean((T_c - T_a)); 

%Heat Gain
Q_i = tau*alpha*sum(I)*A; 

%Total Heat
F = 1; %collecotr efficiency factor
Q_H = F*(Q_i - Q_o);

%Plotting
plot(A, Q_H)
xlabel('Area m^2');
ylabel('Q_H (W)');

