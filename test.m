%% Parameters
m_total = 1095.*10.^12;                                                     %estimated carbon stocks in permafrost
A = integral(@(x) (6371.*1000).^2.*2.*pi.*sin(x), 0, pi/6);                 %Arctic surface area N60 deg (3.4168e+13 m^2)
Z_L = 15;                                                                   %max depth of surface slab
cdensity = @(z) m_total./(A.*Z_L);                                          %density function for carbon storage, uniform

options = optimset('Display','off');                                        %fsolve options
%% Comments
%{ 
 - if the surface temperature is smaller than the constant temperature at Z_L,
the temperature gradient equation makes is such that temperatures increase
as you move deeper.  This implies the depth where temperature is zero is
actually below the bottom of our surface slab. We use an if else block to
account for this 

- the values (9.75, 104, 201.73, 407.4, 0.6, 1875) correspond to modern day
levels for FO, FA, Q, mu, delta and nu -- nu being concentration of CH4

- rate equation for methane does not yet include sink terms 

%}

%% Solve EBM and extend results to Arctic 

x0 = 1; %initial point for fsolve to start at; 1< warm bias, 1> cold bias

%equilibrium temp under modern, globally-averaged forcings 
eq = fsolve(@(tau) eqtemp(tau, 9.75, 104, 201.73, 407.4, 0.6, 1875),x0, options); %equilibrium solution to EBM
Z_0 = fsolve(@(z) tempgrad(eq, 9.75, 104, 201.73, 407.4, 0.6, 1875, z)-1,0, options); %depth where temperature is zero 

%compute mass per unit area 
m_lossCO2_1 = integral(@(z) rateCO2(tempgrad(eq, 9.75, 104, 201.73, 407.4, 0.6, 1875, z)).*cdensity(z), Z_0, 0, 'ArrayValued', true); %carbon released in unit area as CO2
m_lossCH4_1 = integral(@(z) rateCH4(tempgrad(eq, 9.75, 104, 201.73, 407.4, 0.6, 1875, z)).*cdensity(z), Z_0, 0, 'ArrayValued', true); %carbon released in unit area as CH4


m_totalCO2_1 = m_lossCO2_1.*A %mass per unit area times total arctic area 
m_totalCH4_1 = m_lossCH4_1.*A %mass per unit area times total arctic area 

new_concentrations_1 = concentration([m_totalCO2_1 m_totalCH4_1], [407.4 1875])
%% Integrate over latitudes method
%define latitude dependent insolation function from McGehee
relinsol = @(x,y) 2.*sqrt(1-(sqrt(1-sin(y).^2).*sin(0.4101524).*cos(x)-sin(y).*cos(0.4101524)).^2)/pi.^2;
coeff = @(y) integral(@(x) relinsol(x,y), 0, 2.*pi);
insol = @(y) 343.*coeff(y);

%function to determine equilibrium temp and zero temp depth
eqf = @(y) fsolve(@(tau) eqtemp(tau, 9.75, 104, insol(y), 407.4, 0.6, 1875),x0, options);
Z = @(y) fsolve(@(z) tempgrad(eqf(y), 9.75, 104, insol(y), 407.4, 0.6, 1875, z)-1,0, options);

%functions for thaw mass at a latitude
mCO2 = @(y) integral(@(z) rateCO2(tempgrad(eqf(y), 9.75, 104, insol(y), 407.4, 0.6, 1875, z)).*cdensity(z), Z(y), 0,'ArrayValued', true);
mCH4 = @(y) integral(@(z) rateCH4(tempgrad(eqf(y), 9.75, 104, insol(y), 407.4, 0.6, 1875, z)).*cdensity(z), Z(y), 0,'ArrayValued', true);

%integrate over latitudes
m_totalCO2_2 = integral(@(y) 2.*pi.*(6371.*1000).^2.*cos(y).*mCO2(y), pi/3, pi/2)
m_totalCH4_2 = integral(@(y) 2.*pi.*(6371.*1000).^2.*cos(y).*mCH4(y), pi/3, pi/2)

new_concentrations_2 = concentration([m_totalCO2_2 m_totalCH4_2], [407.4 1875])
%%

function r2 = Q10(r1, t1, t2, sensitivity)
    r2 = r1.*(sensitivity.^((t2-t1)./10));
end
function r = rateCO2(tau)
%{
Function for scaling rate of decomposition according to Q10 sensitivity
Args: tau - nondimensional temperature 
Returns: r - decomp rate CO2

current function uses Alternative #2 methodology 
%}

%Constants
gammaA_ms = 0.013; %fraction of carbon in active pool, mineral 
gammaS_ms = 0.107; %fraction of carbon in slow pool, mineral 
gammaA_o = 0.015; %fraction of carbon in active pool, organic 
gammaS_o = 0.293; %fraction of carbon in slow pool, organic
fms = 0.708; %mass fraction of carbon in mineral soils
fo = 0.292; %mass fraction of carbon in organic soils  
A_msan = 0.05; %Area fraction in mineral with anaerobic decomp
A_oan = 0.8; %Area fraction in organic with anaerobic decomp
chi_ms = 0.25; %oxidation fraction mineral under anaerobic
chi_o = 0.6; %oxidation fraction organic under anaerobic
R_ana = 0.1; %ratio of decomp speeds anaerobic:aerobic 
C_tot = 1035*10.^12; %total carbon stocks in permafrost 
Q10_a = 1.5; %Q10 sensitivity aerobic
Q10_an = 3; %Q10 sensitivity anaerobic 
k_a_ms = 1./((0.48+0.21)./2); %decomp rate active, mineral soils
k_s_ms = 1./((8.76+6.42)./2); %decomp rate slow, mineral soils
k_a_o = 1./0.41; %decom rate active, organic soils 
k_s_o = 1./7.21; %decomp rate slow, organic soils 

%define rate scaling functions
T = tau.*273.15-273.15;
ratescale_a_ms = @(T) Q10(k_a_ms, 5, T, Q10_a);
ratescale_s_ms = @(T) Q10(k_s_ms, 5, T, Q10_a);
ratescale_a_o = @(T) Q10(k_a_o, 5, T, Q10_a);
ratescale_s_o = @(T) Q10(k_s_o, 5, T, Q10_a);

if tau>1
    rate_co2 = (gammaA_ms.*ratescale_a_ms(T) + gammaS_ms.*ratescale_s_ms(T)).*fms.*((1-A_msan)+ R_ana.*A_msan.*chi_ms)+(gammaA_o.*ratescale_a_o(T) + gammaS_o.*ratescale_s_o(T)).*fo.*((1-A_oan)+ R_ana.*A_oan.*chi_o);
else 
    rate_co2 = 0;
end

r = rate_co2;
end
function r = rateCH4(tau)
%{
Function for scaling rate of decomposition according to Q10 sensitivity
Args: tau - nondimensional temperature 
Returns: r - decomp rate CH4

current function uses Alternative #2 methodology 
%}

%Constants
gammaA_ms = 0.013; %fraction of carbon in active pool, mineral 
gammaS_ms = 0.107; %fraction of carbon in slow pool, mineral 
gammaA_o = 0.015; %fraction of carbon in active pool, organic 
gammaS_o = 0.293; %fraction of carbon in slow pool, organic
fms = 0.708; %mass fraction of carbon in mineral soils
fo = 0.292; %mass fraction of carbon in organic soils  
A_msan = 0.05; %Area fraction in mineral with anaerobic decomp
A_oan = 0.8; %Area fraction in organic with anaerobic decomp
chi_ms = 0.25; %oxidation fraction mineral under anaerobic
chi_o = 0.6; %oxidation fraction organic under anaerobic
R_ana = 0.1; %ratio of decomp speeds anaerobic:aerobic 
C_tot = 1035*10.^12; %total carbon stocks in permafrost 
Q10_a = 1.5; %Q10 sensitivity aerobic
Q10_an = 3; %Q10 sensitivity anaerobic 
k_a_ms = 1./((0.48+0.21)./2); %decomp rate active, mineral soils
k_s_ms = 1./((8.76+6.42)./2); %decomp rate slow, mineral soils
k_a_o = 1./0.41; %decom rate active, organic soils 
k_s_o = 1./7.21; %decomp rate slow, organic soils 

%define rate scaling functions
T = tau.*273.15-273.15;
ratescale_a_ms = @(T) Q10(k_a_ms, 5, T, Q10_a);
ratescale_s_ms = @(T) Q10(k_s_ms, 5, T, Q10_a);
ratescale_a_o = @(T) Q10(k_a_o, 5, T, Q10_a);
ratescale_s_o = @(T) Q10(k_s_o, 5, T, Q10_a);

if tau>1
    rate_ch4 = (gammaA_ms.*ratescale_a_ms(T) + gammaS_ms.*ratescale_s_ms(T)).*R_ana.*(fms.*A_msan.*(1-chi_ms)) + (gammaA_o.*ratescale_a_o(T) + gammaS_o.*ratescale_s_o(T)).*R_ana.*(fo.*A_oan.*(1-chi_o));
else 
    rate_ch4 = 0;
end 

r = rate_ch4;
end

%{
tempgrad function recomputes everything in eqtemp but using the
equilibrium temperature.  Both functions can be combined but its more
legible for the time being to keep them apart. These functions are based on
the heat equation solution with surface flux as a boundary condition 
%}
function F = eqtemp(tau, F_O, F_A, Q, mu, delta, nu)
% modified EBM

% CONSTANTS
T_R = 273.15;                                   % K                     % reference temperature, equal to 0 Celcius
betaconst = 0.63;                               % nondimensional        % fraction of atmospheric irradiance hitting surface           
sigma = 5.670367*10^-8;                         % W/(m^2K^4)            % Stefan-Boltzmann constant 
omega = 0.01;                                   % nondimensional        % albedo switch function smoothness
Gamma = 6.49*10^-3;                             % K/m                   % standard ICAO lapse rate for the troposphere
gamma = Gamma/T_R;                              % 1/m                   % scaled lapse rate for non-dimensional system
k_C = 0.0694;                                   % m^2/kg                % absorption cross-section per unit mass of atmosphere for carbon dioxide
k_W = 0.05905;                                  % m^2/kg                % absorption cross-section per unit mass of atmosphere for water vapour
k_ch4 = 1.4022;                                 % m^2/kg                % absorption cross-section per unit mass of atmosphere for methane 
Z_P = 9000;                                     % m                     % troposphere height      
L_v = 2.2558*10^6;                              % m^2/s^2               % latent heat of vapourization for water
R = 8.3144598;                                  % J/(mol K)             % universal gas constant
R_W = 461.4;                                    % m^2/(s^2 K)           % specific gas constant of water vapour
alpha_w = 0.08;                                 % nondimensional        % warm (ocean) albedo
alpha_c = 0.7;                                  % nondimensional        % cold (ice) albedo                                 % W/m^2                 % average incoming solar radiation above 60N
q = Q/(sigma*T_R^4);                            % nondimensional        % scaled incoming solar radiation at 70N
c_P = 1004;                                     % J/(K kg)              % specific heat capacity of dry atmopshere at T_R
C_LH = 1.0*10^-3;                               % nondimensional        % latent heat bulk exchange coefficient 
C_SH = 1.0*10^-3;                               % nondimensional        % sensible heat bulk exchange coefficient
U_r = 4;                                        % m/s                   % average horizontal wind speed
z_r = 600;                                      % m                     % average height of planetary boundary layer
P_sat = 611.2;                                  % Pa = kg/(m s^2)       % saturated vapour pressure at T_R
rho_wsat = 4.849*10^-3;                         % kg/m^3                % saturated vapour density at T_R
P_0 = 101325.00;                                % Pa                    % pressure at surface 
M =  0.0289644;                                 % kg/mol                % molecular weight of dry air
G_W1 = L_v/(R_W*T_R);                           % nondimensional        % atmospheric water vapour absorption coefficient 1

%%%% new terms for permafrost %%%%
k_L = 2;                                        % W m^-1 K^-1           % thermal conductivity 
Z_L = 15;                                       % m                     % surface slab depth 
H2 = k_L./(Z_L.*sigma.*T_R^3);                  % nondimensional 
T_ZL = -2.35;                                   % C                     % ref temp at bottom of slab    
tau_ZL = (T_ZL+273.17)./273.15;                 % nondimensional        % nondimensional ref temp at bottom of slab

% FUNCTIONS
T_S = tau*T_R;                                  % K                     % un-scaled surface temperature for use in SH and LH
% density of atmosphere (kg/m^3) at the top of PBL as a function of surface temperature
rho = P_0.*M./(R.*(T_S-Gamma.*z_r));                        
% vertical sensible heat transport (W/m^2) due to turbulent flux in the PBL
SH = c_P.*rho.*C_SH.*U_r.*(Gamma.*z_r);
% vertical latent heat transport (W/m^2) due to turbulent flux in the PBL
LE = L_v.*C_LH.*U_r./(R_W.*(T_S-Gamma.*z_r)).*P_sat.*(exp(G_W1.*(T_S - T_R)./T_S) - delta.*exp(G_W1.*(T_S - Gamma.*z_r -T_R)./(T_S - Gamma.*z_r)));
% combined vertical turbulent heat transport in the PBL
F_C = LE + SH;

f_A  = F_A./(sigma.*T_R.^4);                    % nondimensional        % scaled horizontal atmospheric heat transport
f_O  = F_O./(sigma.*T_R.^4);                    % nondimensional        % scaled horizontal aoceanic heat transport
f_C = F_C./(sigma.*T_R.^4);                     % nondimensional        % scaled vertical turbulent heat flux heat transport

% albedo switch function (nondimensional)
a = 1./2.*((alpha_w + alpha_c)+(alpha_w-alpha_c).*tanh((tau-1)./omega));    

%modified eta_cl
eta_Cl = 0.3736;                                % nondimensional        % absorption due to clouds
G_c = 1.52./(10.^6).*k_C.*1.03.*10.^4;          % mol/micromol          % atmospheric carbon dioxide absorption coefficient
eta_C1 = 1 - exp(-G_c.*mu);                     % nondimensional        % absorption due to carbon dioxide         
G_ch = 0.554./(10.^9).*k_ch4.*1.03.*10.^4;
eta_ch = 1 - exp(-G_ch.*nu);
G_W2 = k_W.*rho_wsat./gamma;                    % nondimensional        % atmospheric water vapour absorption coefficient 2
wvinteg1 = @(w) 1./w.*exp(G_W1.*(w-1)./w);      % nondimensional        % temperature-dependent part of Clausius-Clapeyron equation integrated over troposphere
eta_W1 = 1 - exp(-delta.*G_W2.*integral(wvinteg1, tau-(gamma.*Z_P), tau));   % nondimensional        % absorption due to water vapour     

eta = 1 - (1-eta_C1).*(1-eta_W1).*(1-eta_Cl)*(1-eta_ch);   % nondimensional        % total atmospheric longwave absorption

F = tau_ZL - tau + (1).*(f_O - (1-betaconst).*f_C + (1-a).*q.*(1 - 0.2324 - 0.1212) - (1-betaconst.*eta).*tau.^4 + betaconst.*(f_A + 0.2324.*q))./H2;
end
function T = tempgrad(eqtau, F_O, F_A, Q, mu, delta, nu, z)
%describes temperature as a function of depth, z (m), at equilibrium

% CONSTANTS
T_R = 273.15;                                   % K                     % reference temperature, equal to 0 Celcius
betaconst = 0.63;                               % nondimensional        % fraction of atmospheric irradiance hitting surface           
sigma = 5.670367*10^-8;                         % W/(m^2K^4)            % Stefan-Boltzmann constant 
omega = 0.01;                                   % nondimensional        % albedo switch function smoothness
Gamma = 6.49*10^-3;                             % K/m                   % standard ICAO lapse rate for the troposphere
gamma = Gamma/T_R;                              % 1/m                   % scaled lapse rate for non-dimensional system
k_C = 0.0694;                                  % m^2/kg                % absorption cross-section per unit mass of atmosphere for carbon dioxide
k_W = 0.05905;                                  % m^2/kg                % absorption cross-section per unit mass of atmosphere for water vapour
k_ch = 1.4022;
Z_P = 9000;                                     % m                     % troposphere height      
L_v = 2.2558*10^6;                              % m^2/s^2               % latent heat of vapourization for water
R = 8.3144598;                                  % J/(mol K)             % universal gas constant
R_W = 461.4;                                    % m^2/(s^2 K)           % specific gas constant of water vapour
alpha_w = 0.08;                                 % nondimensional        % warm (ocean) albedo
alpha_c = 0.7;                                  % nondimensional        % cold (ice) albedo                                 % W/m^2                 % average incoming solar radiation above 60N
q = Q/(sigma*T_R^4);                            % nondimensional        % scaled incoming solar radiation at 70N
c_P = 1004;                                     % J/(K kg)              % specific heat capacity of dry atmopshere at T_R
C_LH = 1.0*10^-3;                               % nondimensional        % latent heat bulk exchange coefficient 
C_SH = 1.0*10^-3;                               % nondimensional        % sensible heat bulk exchange coefficient
U_r = 4;                                        % m/s                   % average horizontal wind speed
z_r = 600;                                      % m                     % average height of planetary boundary layer
P_sat = 611.2;                                  % Pa = kg/(m s^2)       % saturated vapour pressure at T_R
rho_wsat = 4.849*10^-3;                         % kg/m^3                % saturated vapour density at T_R
P_0 = 101325.00;                                % Pa                    % pressure at surface 
M =  0.0289644;                                 % kg/mol                % molecular weight of dry air
G_W1 = L_v/(R_W*T_R);                           % nondimensional        % atmospheric water vapour absorption coefficient 1

%%%% new terms %%%%
k_L = 2;
Z_L = 15;
zeta = z./Z_L;
H2 = k_L./(Z_L.*sigma.*T_R^3);
T_ZL = -2.35;
tau_ZL = (T_ZL+273.17)./273.15;

% FUNCTIONS
T_S = eqtau*T_R;                                  % K                     % un-scaled surface temperature for use in SH and LH
% density of atmosphere (kg/m^3) at the top of PBL as a function of surface temperature
rho = P_0.*M./(R.*(T_S-Gamma.*z_r));                        
% vertical sensible heat transport (W/m^2) due to turbulent flux in the PBL
SH = c_P.*rho.*C_SH.*U_r.*(Gamma.*z_r);
% vertical latent heat transport (W/m^2) due to turbulent flux in the PBL
LE = L_v.*C_LH.*U_r./(R_W.*(T_S-Gamma.*z_r)).*P_sat.*(exp(G_W1.*(T_S - T_R)./T_S) - delta.*exp(G_W1.*(T_S - Gamma.*z_r -T_R)./(T_S - Gamma.*z_r)));
% combined vertical turbulent heat transport in the PBL
F_C = LE + SH;

f_A  = F_A./(sigma.*T_R.^4);                    % nondimensional        % scaled horizontal atmospheric heat transport
f_O  = F_O./(sigma.*T_R.^4);                    % nondimensional        % scaled horizontal aoceanic heat transport
f_C = F_C./(sigma.*T_R.^4);                     % nondimensional        % scaled vertical turbulent heat flux heat transport

% albedo switch function (nondimensional)
a = 1./2.*((alpha_w + alpha_c)+(alpha_w-alpha_c).*tanh((eqtau-1)./omega));    

%modified eta_cl
eta_Cl = 0.3736;                                % nondimensional        % absorption due to clouds
G_c = 1.52./(10.^6).*k_C.*1.03.*10.^4;          % mol/micromol          % atmospheric carbon dioxide absorption coefficient
eta_C1 = 1 - exp(-G_c.*mu);                     % nondimensional        % absorption due to carbon dioxide         
G_ch = 0.554./(10.^9).*k_ch.*1.03.*10.^4;
eta_ch = 1 - exp(-G_ch.*nu);
G_W2 = k_W.*rho_wsat./gamma;                    % nondimensional        % atmospheric water vapour absorption coefficient 2
wvinteg1 = @(w) 1./w.*exp(G_W1.*(w-1)./w);      % nondimensional        % temperature-dependent part of Clausius-Clapeyron equation integrated over troposphere
eta_W1 = 1 - exp(-delta.*G_W2.*integral(wvinteg1, eqtau-(gamma.*Z_P), eqtau));   % nondimensional        % absorption due to water vapour     

eta = 1 - (1-eta_C1).*(1-eta_W1).*(1-eta_Cl)*(1-eta_ch);   % nondimensional        % total atmospheric longwave absorption

T = tau_ZL + (zeta+1).*(f_O - (1-betaconst).*f_C + (1-a).*q.*(1 - 0.2324 - 0.1212) - (1-betaconst.*eta).*eqtau.^4 + betaconst.*(f_A + 0.2324.*q))./H2;
end
function c = concentration(masses, initial_concentrations)
%{
Computes new concentrations for CO2 and CH4 given the mass of each gas 
input to the atmosphere and the previous global concentrations.

Args:  masses - tuple containing added mass of carbon in the form of CO2 
                and CH4 ([mCO2 mCH4])
       initial_concentrations - tuple of concentrations before permafrost
                                respiration, [mu_i nu_i] in units [ppm ppb]

Returns: c - tuple of new global atmospheric concentrations, 
             [mu_f nu_f] in units [ppm ppb]
%}

cco2_i = initial_concentrations(1); %initial concentration CO2 (ppm)
cch4_i = initial_concentrations(2); %initial concentration CH4 (ppb)

MCO2 = 44.01;   %Molar masses of CO2, CH4, C, O, H and atmosphere 
MCH4 = 16.04;
MC = 12.0107;
MO = 15.999;
MH = 1.00784;
M = 28.97; 

mco2C = masses(1); %mass of carbon respired by permafrost as CO2
mch4C = masses(2); %mass of carbon respired by permafrost as CH4
mco2O = (mco2C./MC).*2.*MO; %take into account mass of oxygen needed in formation of CO2
mch4H = (mch4C./MC).*4.*MH; %take into account mass of hydrogen needed in formation of CH4
mco2added = mco2C + mco2O;  %total mass of CO2 respired 
mch4added = mch4C + mch4H;  %total mass of CH4 respired 

cco2 = cco2_i./(10.^6);  %ppm to fraction
cch4 = cch4_i./(10.^9);  %ppb to fraction
rho = 1.2;   %atmospheric density 

%atmosphere volume calculation
radius_troposphere = @(x) 16000.*x./pi + 9000 + 6371.*1000;
volume = 2.*integral2(@(x,y) (1/3).*sin(x).*(radius_troposphere(x).^3 - (6371.*1000).^3), 0, pi./2, 0, 2.*pi);

M_atm = rho.*volume; %total mass of atmosphere
mco2initial = (MCO2./M).*cco2.*M_atm; %initial mass of CO2 in atmosphere using mol to mass ratio
mch4initial = (MCH4./M).*cch4.*M_atm; %initial mass of CH4 in atmosphere

M_atmnew = M_atm + mco2added + mch4added;   %new mass of atmosphere
mco2new = mco2initial + mco2added;  %new mass of CO2 in atmosphere
mch4new = mch4initial + mch4added;  %new mass of CH4 in atmosphere

co2ratio = (M./MCO2).*(mco2new./M_atmnew);  %mass concentration to mol concentration
ch4ratio = (M./MCH4).*(mch4new./M_atmnew);  

c(1) = co2ratio.*10.^6; %fraction to ppm
c(2) = ch4ratio.*10.^9; %fraction to ppb
end