fsolve(@(t) equilibrium(t, 9.75, 104, 650, 0.6), 0.9)
%equilibrium(1.018, 9.75, 104, 650, 0.6)

function F = equilibrium(tau, F_O, F_A, mu, delta)
T_R = 273.15;                                   % K                     % reference temperature, equal to 0 Celcius
betaconst = 0.63;                               % nondimensional        % fraction of atmospheric irradiance hitting surface           
sigma = 5.670367*10^-8;                         % W/(m^2K^4)            % Stefan-Boltzmann constant 
omega = 0.01;                                   % nondimensional        % albedo switch function smoothness
Gamma = 6.49*10^-3;                             % K/m                   % standard ICAO lapse rate for the troposphere
gamma = Gamma/T_R;                              % 1/m                   % scaled lapse rate for non-dimensional system
k_C = 0.07424;                                  % m^2/kg                % absorption cross-section per unit mass of atmosphere for carbon dioxide
k_W = 0.05905;                                  % m^2/kg                % absorption cross-section per unit mass of atmosphere for water vapour
Z_P = 9000;                                     % m                     % troposphere height      
L_v = 2.2558*10^6;                              % m^2/s^2               % latent heat of vapourization for water
R = 8.3144598;                                  % J/(mol K)             % universal gas constant
R_W = 461.4;                                    % m^2/(s^2 K)           % specific gas constant of water vapour
alpha_w = 0.08;                                 % nondimensional        % warm (ocean) albedo
alpha_c = 0.7;                                  % nondimensional        % cold (ice) albedo
Q = 201.73;                                      % W/m^2                 % incoming solar radiation at 70N
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
rho_O = 1000;
c_O = 4281;
Z_O = 70;
c_P = 1004;
rho_L = 1224;
c_L = 2050;
Z_L = 15;
K_L = 2;



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

eta_Cl = 0.3729;                                % nondimensional        % absorption due to clouds
G_c = 1.52./(10.^6).*k_C.*1.03.*10.^4;          % mol/micromol          % atmospheric carbon dioxide absorption coefficient
eta_C1 = 1 - exp(-G_c.*mu);                     % nondimensional        % absorption due to carbon dioxide         
G_W2 = k_W.*rho_wsat./gamma;                    % nondimensional        % atmospheric water vapour absorption coefficient 2
wvinteg1 = @(w) 1./w.*exp(G_W1.*(w-1)./w);      % nondimensional        % temperature-dependent part of Clausius-Clapeyron equation integrated over troposphere
eta_W1 = 1 - exp(-delta.*G_W2.*integral(wvinteg1, tau-(gamma.*Z_P), tau));   % nondimensional        % absorption due to water vapour     
eta = 1 - (1-eta_C1).*(1-eta_W1).*(1-eta_Cl);   % nondimensional        % total atmospheric longwave absorption

%temperature gradient, surface flux and tau_bottom boundary conditions 
H2 = K_L./(Z_L.*sigma.*T_R.^3);
tau_bottom = (273.15-2.35)./273.15;

zeta0 = 1 - (1-tau_bottom).*H2./(f_O - (1-betaconst).*f_C + (1-a).*q.*(1 - 0.2324 - 0.1212) - (1-betaconst.*eta).*tau.^4 + betaconst.*(f_A+0.2324.*q))
net_flux = f_O - (1-betaconst).*f_C + (1-a).*q.*(1 - 0.2324 - 0.1212) - (1-betaconst.*eta).*tau.^4 + betaconst.*(f_A+0.2324.*q)
F = tau_bottom - tau + ((1-0)./H2).*(f_O - (1-betaconst).*f_C + (1-a).*q.*(1 - 0.2324 - 0.1212) - (1-betaconst.*eta).*tau.^4 + betaconst.*(f_A+0.2324.*q));

end 

