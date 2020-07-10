classdef constants
    properties (Constant)
        % ----- ORIGINAL EBM CONSTANTS ---- %
        T_R = 273.15;                                   % K                     % reference temperature, equal to 0 Celcius
        betaconst = 0.63;                               % nondimensional        % fraction of atmospheric irradiance hitting surface           
        sigma = 5.670367*10^-8;                         % W/(m^2K^4)            % Stefan-Boltzmann constant 
        omega = 0.01;                                   % nondimensional        % albedo switch function smoothness
        Gamma = 6.49*10^-3;                             % K/m                   % standard ICAO lapse rate for the troposphere
        gamma = constants.Gamma/constants.T_R;          % 1/m                   % scaled lapse rate for non-dimensional system
        k_C = 0.0694;                                   % m^2/kg                % absorption cross-section per unit mass of atmosphere for carbon dioxide
        k_W = 0.05905;                                  % m^2/kg                % absorption cross-section per unit mass of atmosphere for water vapour
        k_M = 1.4022;                                   % m^2/kg                % absorption cross-section per unit mass of atmosphere for methane 
        Z_P = 9000;                                     % m                     % troposphere height      
        L_v = 2.2558*10^6;                              % m^2/s^2               % latent heat of vapourization for water
        R = 8.3144598;                                  % J/(mol K)             % universal gas constant
        R_W = 461.4;                                    % m^2/(s^2 K)           % specific gas constant of water vapour
        alpha_w = 0.08;                                 % nondimensional        % warm (ocean) albedo
        alpha_c = 0.7;                                  % nondimensional        % cold (ice) albedo                                 % W/m^2                
        c_P = 1004;                                     % J/(K kg)              % specific heat capacity of dry atmopshere at T_R
        C_LH = 1.0*10^-3;                               % nondimensional        % latent heat bulk exchange coefficient 
        C_SH = 1.0*10^-3;                               % nondimensional        % sensible heat bulk exchange coefficient
        U_r = 4;                                        % m/s                   % average horizontal wind speed
        z_r = 600;                                      % m                     % average height of planetary boundary layer
        P_sat = 611.2;                                  % Pa = kg/(m s^2)       % saturated vapour pressure at T_R
        rho_wsat = 4.849*10^-3;                         % kg/m^3                % saturated vapour density at T_R
        P_0 = 101325.00;                                % Pa                    % pressure at surface 
        M =  0.0289644;                                 % kg/mol                % molecular weight of dry air
        G_W1 = constants.L_v/(constants.R_W*constants.T_R); % nondimensional        % atmospheric water vapour absorption coefficient 1

        %----- PERMAFROST CONSTANTS -----%
        A_tot = integral(@(x) (6371.*1000).^2.*2.*pi.*sin(x), 0, pi/6);
        A_p = 16.*10.^6.*1000.^2;
        Ctot = 1.095.*10.^15;                           % kg                    %estimated carbon stocks
        
        k_L = 2;                                        % W m^-1 K^-1           % thermal conductivity 
        Z_L = 15;                                       % m                     % surface slab depth 
        H2 = constants.k_L./(constants.Z_L.*constants.sigma.*constants.T_R^3); % nondimensional 
        T_ZL = -10.35;                                   % C                     % ref temp at bottom of slab    
        tau_ZL = (constants.T_ZL+constants.T_R)./constants.T_R; % nondimensional        % nondimensional ref temp at bottom of slab
        
    
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
        Q10_a = 1.5; %Q10 sensitivity aerobic
        Q10_an = 3; %Q10 sensitivity anaerobic 
        k_a_ms = 1./((0.48+0.21)./2); %decomp rate active, mineral soils
        k_s_ms = 1./((8.76+6.42)./2); %decomp rate slow, mineral soils
        k_a_o = 1./0.41; %decom rate active, organic soils 
        k_s_o = 1./7.21; %decomp rate slow, organic soils 
        
        decay_CH4 = log(2)./9.1;
        
        MCO2 = 44.01e-3;   %Molar masses of CO2, CH4 and C kg mol^-1
        MCH4 = 16.04e-3;
        MC = 12.0107e-3;
    
    end 
end 