%should output derivative of our state variable and mass should be updated 
%in loop during the AB3 process

classdef ebm < dynamicprops
    properties  %DECLARE PROPERTIES OF OUR EBM OBJECT 
        %forcings 
        F_O {mustBeNumeric}
        F_A {mustBeNumeric}
        Q {mustBeNumeric}
        mu {mustBeNumeric}
        nu {mustBeNumeric}
        delta {mustBeNumeric}
        
        %properties
        tau_s {mustBeNumeric}
        state %(mass dist, co2cum_emissions, ch4cum_emissions)
         
    end 
    
    methods
        %init
        function e = ebm(F_O, F_A, Q, delta)
            import constants.*;
            
            %forcings
            e.F_O = F_O;
            e.F_A = F_A;
            e.Q = Q;
            e.delta = delta;
            
            %default params (so you can test functions right away)
            e.mu = 407.4;
            e.nu = 1875;
            
            %initial carbon distribution, 30 pools uniform 
            e.addprop('P')
            e.P = 30;
            Cdist0 = (constants.Ctot./constants.A_p).*ones([1 e.P])./e.P;
            
            %initial ghg emissions 
            CO20 = 0;
            CH40 = 0;
            emissions0 = [CO20 CH40];
            
            %initial state variable 
            e.state = horzcat(Cdist0, emissions0);
        end 
        %equilibrium surface temperature
        function eqtemp(e,x0)
            %calculate atmospheric concentrations 
            conc_new = ebm.concentration([e.state(end-1) e.state(end)], [e.mu e.nu]);
            e.mu = conc_new(1);
            e.nu = conc_new(2);
            
            %solve implicit equation to find tau_s
            e.tau_s = fsolve(@(t) ebm.eq_implicit(t, e.F_O, e.F_A, e.Q, e.mu, e.nu, e.delta),x0);
        end 
        %equilibrium temperature profile
        function T = eqtempProfile(e,z)
            T = ebm.tempgrad(e.tau_s, e.F_O, e.F_A, e.Q, e.mu, e.nu, e.delta, z);
        end 
        %computes emissions and mass loss - d(state)/dt
        function d = dynamics(e)
            %figure out where zero temp depth is 
            zero_depth = fsolve(@(z) eqtempProfile(e,z)-1,0);
            
            if zero_depth>-constants.Z_L & zero_depth<0
                %figure out which index to increment to 
                delta_z = -constants.Z_L./e.P;
                zero_index = round(zero_depth./delta_z);

                %iterate through depths, compute carbon emissions as CH4 and CO2
                %reimann sum 
                for i = 1:zero_index
                    c(i) = e.state(i).*ebm.HRCO2(eqtempProfile(e,-(i-0.5)*constants.Z_L./e.P), (i-0.5)./e.P);
                    m(i) = e.state(i).*ebm.HRCH4(eqtempProfile(e,-(i-0.5)*constants.Z_L./e.P), (i-0.5)./e.P);
                    losses(i) = -(c(i) + m(i));    %net carbon losses in a layer
                end 
                
                %account for depths with no decay for state derivative
                missing = zeros([1, e.P-zero_index]);
                losses = horzcat(losses, missing);

                %now account for CO2 and CH4 specifics
                mCO2 = (constants.MCO2./constants.MC).*constants.A_p.*sum(c) + (constants.MCO2./constants.MCH4).*constants.decay_CH4.*e.state(end);
                mCH4 = (constants.MCH4./constants.MC).*constants.A_p.*sum(m) - constants.decay_CH4.*e.state(end);
            
                emissions = [mCO2 mCH4];
                d = horzcat(losses, emissions);
            else 
                d = zeros([1 32]);
            end 
        end  
    end 
    
    methods (Static)    %CHANGE ACCESS TO PRIVATE LATER
        %generic functions
        function F = eq_implicit(tau, F_O, F_A, Q, mu, nu, delta)
            T_S = tau*constants.T_R;                             
            % density of atmosphere (kg/m^3) at the top of PBL as a function of surface temperature
            rho = constants.P_0.*constants.M./(constants.R.*(T_S-constants.Gamma.*constants.z_r));                        
            % vertical sensible heat transport (W/m^2) due to turbulent flux in the PBL
            SH = constants.c_P.*rho.*constants.C_SH.*constants.U_r.*(constants.Gamma.*constants.z_r);
            % vertical latent heat transport (W/m^2) due to turbulent flux in the PBL
            LE = constants.L_v.*constants.C_LH.*constants.U_r./(constants.R_W.*(T_S-constants.Gamma.*constants.z_r)).*constants.P_sat.*(exp(constants.G_W1.*(T_S - constants.T_R)./T_S) - delta.*exp(constants.G_W1.*(T_S - constants.Gamma.*constants.z_r - constants.T_R)./(T_S - constants.Gamma.*constants.z_r)));
            % combined vertical turbulent heat transport in the PBL
            F_C = LE + SH;

            f_A  = F_A./(constants.sigma.*constants.T_R.^4);                    % nondimensional        % scaled horizontal atmospheric heat transport
            f_O  = F_O./(constants.sigma.*constants.T_R.^4);                    % nondimensional        % scaled horizontal aoceanic heat transport
            f_C = F_C./(constants.sigma.*constants.T_R.^4);                     % nondimensional        % scaled vertical turbulent heat flux heat transport
            q = Q./(constants.sigma.*constants.T_R^4);
            
            % albedo switch function (nondimensional)
            a = 1./2.*((constants.alpha_w + constants.alpha_c)+(constants.alpha_w-constants.alpha_c).*tanh((tau-1)./constants.omega));    

            %eta
            eta_Cl = 0.3736;                               
            G_c = 1.52./(10.^6).*constants.k_C.*1.03.*10.^4;         
            eta_C1 = 1 - exp(-G_c.*mu);                           
            G_M = 0.554./(10.^9).*constants.k_M.*1.03.*10.^4;
            eta_M = 1 - exp(-G_M.*nu);
            G_W2 = constants.k_W.*constants.rho_wsat./constants.gamma;                   
            wvinteg1 = @(w) 1./w.*exp(constants.G_W1.*(w-1)./w);     
            eta_W1 = 1 - exp(-delta.*G_W2.*integral(wvinteg1, tau-(constants.gamma.*constants.Z_P), tau));   
            eta = 1 - (1-eta_C1).*(1-eta_W1).*(1-eta_Cl).*(1-eta_M);   

            %final step 
            F = constants.tau_ZL - tau + (1).*(f_O - (1-constants.betaconst).*f_C + (1-a).*q.*(1 - 0.2324 - 0.1212) - (1-constants.betaconst.*eta).*tau.^4 + constants.betaconst.*(f_A + 0.2324.*q))./constants.H2;
        end 
        function T = tempgrad(tau, F_O, F_A, Q, mu, nu, delta, z)
            zeta = z./constants.Z_L;
       
            T_S = tau*constants.T_R;                             
            % density of atmosphere (kg/m^3) at the top of PBL as a function of surface temperature
            rho = constants.P_0.*constants.M./(constants.R.*(T_S-constants.Gamma.*constants.z_r));                        
            % vertical sensible heat transport (W/m^2) due to turbulent flux in the PBL
            SH = constants.c_P.*rho.*constants.C_SH.*constants.U_r.*(constants.Gamma.*constants.z_r);
            % vertical latent heat transport (W/m^2) due to turbulent flux in the PBL
            LE = constants.L_v.*constants.C_LH.*constants.U_r./(constants.R_W.*(T_S-constants.Gamma.*constants.z_r)).*constants.P_sat.*(exp(constants.G_W1.*(T_S - constants.T_R)./T_S) - delta.*exp(constants.G_W1.*(T_S - constants.Gamma.*constants.z_r - constants.T_R)./(T_S - constants.Gamma.*constants.z_r)));
            % combined vertical turbulent heat transport in the PBL
            F_C = LE + SH;

            f_A  = F_A./(constants.sigma.*constants.T_R.^4) ;                   % nondimensional        % scaled horizontal atmospheric heat transport
            f_O  = F_O./(constants.sigma.*constants.T_R.^4);                    % nondimensional        % scaled horizontal aoceanic heat transport
            f_C = F_C./(constants.sigma.*constants.T_R.^4);                     % nondimensional        % scaled vertical turbulent heat flux heat transport
            q = Q./(constants.sigma.*constants.T_R^4);
            
            % albedo switch function (nondimensional)
            a = 1./2.*((constants.alpha_w + constants.alpha_c)+(constants.alpha_w-constants.alpha_c).*tanh((tau-1)./constants.omega));    

            %eta
            eta_Cl = 0.3736;                               
            G_c = 1.52./(10.^6).*constants.k_C.*1.03.*10.^4;         
            eta_C1 = 1 - exp(-G_c.*mu);                           
            G_M = 0.554./(10.^9).*constants.k_M.*1.03.*10.^4;
            eta_M = 1 - exp(-G_M.*nu);
            G_W2 = constants.k_W.*constants.rho_wsat./constants.gamma;                   
            wvinteg1 = @(w) 1./w.*exp(constants.G_W1.*(w-1)./w);     
            eta_W1 = 1 - exp(-delta.*G_W2.*integral(wvinteg1, tau-(constants.gamma.*constants.Z_P), tau));   
            eta = 1 - (1-eta_C1).*(1-eta_W1).*(1-eta_Cl).*(1-eta_M);
            
            %Final step
            T = constants.tau_ZL + (zeta+1).*(f_O - (1-constants.betaconst).*f_C + (1-a).*q.*(1 - 0.2324 - 0.1212) - (1-constants.betaconst.*eta).*tau.^4 + constants.betaconst.*(f_A + 0.2324.*q))./constants.H2;
        end 
        function r2 = Q10(r1, t1, t2, sensitivity)
            r2 = r1.*(sensitivity.^((t2-t1)./10));
        end
        function r = HRCO2(tau, zeta)
            T = tau.*constants.T_R-constants.T_R;
            ratescale_a_ms = @(T) ebm.Q10(constants.k_a_ms, 5, T, constants.Q10_a);
            ratescale_s_ms = @(T) ebm.Q10(constants.k_s_ms, 5, T, constants.Q10_a);
            ratescale_a_o = @(T) ebm.Q10(constants.k_a_o, 5, T, constants.Q10_a);
            ratescale_s_o = @(T) ebm.Q10(constants.k_s_o, 5, T, constants.Q10_a);


            rate_co2 = @(t) (constants.gammaA_ms.*ratescale_a_ms(t) + constants.gammaS_ms.*ratescale_s_ms(t)).*constants.fms.*((1-constants.A_msan)+ constants.R_ana.*constants.A_msan.*constants.chi_ms)+(constants.gammaA_o.*ratescale_a_o(t) + constants.gammaS_o.*ratescale_s_o(t)).*constants.fo.*((1-constants.A_oan)+ constants.R_ana.*constants.A_oan.*constants.chi_o);
            %account for seasonality over a year
            months = linspace(0,11,12);
            deltaT = 5;
            for i=1:12
                if T+(1-zeta).*deltaT.*sin(months(i)*2*pi./11) > 0
                    rates(i) = rate_co2(T+(1-zeta).*deltaT.*sin(months(i)*2*pi./11));
                else 
                    rates(i) = 0;
                end
            end 
            r = mean(rates);
        end 
        function r = HRCH4(tau, zeta)
            T = tau.*constants.T_R-constants.T_R;
            ratescale_a_ms = @(T) ebm.Q10(constants.k_a_ms, 5, T, constants.Q10_a);
            ratescale_s_ms = @(T) ebm.Q10(constants.k_s_ms, 5, T, constants.Q10_a);
            ratescale_a_o = @(T) ebm.Q10(constants.k_a_o, 5, T, constants.Q10_a);
            ratescale_s_o = @(T) ebm.Q10(constants.k_s_o, 5, T, constants.Q10_a);

            rate_ch4 = @(t) (constants.gammaA_ms.*ratescale_a_ms(t) + constants.gammaS_ms.*ratescale_s_ms(t)).*constants.R_ana.*constants.fms.*constants.A_msan.*(1-constants.chi_ms) + (constants.gammaA_o.*ratescale_a_o(t) + constants.gammaS_o.*ratescale_s_o(t)).*constants.R_ana.*constants.fo.*constants.A_oan.*(1-constants.chi_o);
            %account for seasonality over a year
            months = linspace(0,11,12);
            deltaT = 5;
            for i=1:12
                if T+(1-zeta).*deltaT.*sin(months(i)*2*pi./11) > 0
                    rates(i) = rate_ch4(T+(1-zeta).*deltaT.*sin(months(i)*2*pi./11));
                else 
                    rates(i) = 0;
                end
            end 
            r = mean(rates);
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

            MCO2 = 44.01e-3;   %Molar masses of CO2, CH4, C, O, H and atmosphere kg mol^-1
            MCH4 = 16.04e-3;
            MC = 12.0107e-3;
            M = 28.97e-3; 

            mco2added = masses(1);
            mch4added = masses(2);

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
    end 
end 
