%{
The following code is an implementation of a temperature distribution 
mechanism as outlined in the "Global significance of permafrost" document. 
We will compare this method to the more simple idea of scaling results in a
unit area to the entire Arctic through multiplication by Arctic surface
area. We will look at GHG concentration changes over a year due to 
permafrost thaw. 
%}

%% Shared parameters
A = integral(@(x) (6371.*1000).^2.*2.*pi.*sin(x), 0, pi/6);     % Arctic surface area N60 deg (3.4168e+13 m^2)
m = 1300.*10^12;                                                % estimated carbon stocks (1300 Pg)
m_unit = m./A;                                                  % mass of carbon in an area unit (kg/m^2) 
t = 1.018;                                                      % non-dimensional temp 5 degrees for comparison's sake
c_init = [407.4 1875];                                          % initial GHG concentrations [CO2 CH4] ([ppm ppb])                   

%% Simple model 
mCO2_simple = rateCO2(t).*mthaw(t, m_unit).*A;                 % mass of carbon respired as CO2
mCH4_simple = rateCH4(t).*mthaw(t, m_unit).*A;                 % mass of carbon respired as CH4
masses1 = [mCO2_simple mCH4_simple]                       % tuple of masses for compatibility with concentration function
conc1 = concentration(masses1, c_init)

%% Temperature Distribution model
masses2 = integralmethod(t, m_unit)
conc2 = concentration(masses2, c_init)

%% related to permafrost 2 WONT WORK HERE
%Integrate over slabs w/ latitude-dependent insolation over Arctic SA

% use insolation coefficient distribution from McGehee
relinsol = @(x,y) 2.*sqrt(1-(sqrt(1-sin(y).^2).*sin(0.4101524).*cos(x)-sin(y).*cos(0.4101524)).^2)/pi.^2;
coeff = @(y) integral(@(x) relinsol(x,y), 0, 2.*pi);
insol = @(y) 343.*coeff(y);                                                 %average insolation times relative factor, W/m^2

% create functions which are compatible with Arctic SA integral
equilibrium = @(y) fsolve(@(tau) eqtemp(tau, 9.75, 104, insol(y), 407.4, 0.6, 1875),1.03, options); %latitude-dependent eq temp
zero_depth = @(y) fsolve(@(z) tempgrad(equilibrium(y), 9.75, 104, insol(y), 407.4, 0.6, 1875, z)-1,0, options); %zero temp depth as a function of latitude 
mass_carbon = @(y) integral(@(z) rate(tempgrad(equilibrium(y), 9.75, 104, insol(y), 407.4, 0.6, 1875, z)).*cdensity(z), zero_depth(y), 0);

% integrate over Arctic SA
m_totalV2 = integral(@(y) mass_carbon(y).*2.*pi.*(6371.*1000).^2.*cos(y).*cos(y), pi./3, pi./2, 'ArrayValued', true)

%% Summary 
%{
Under the assumption that the mean temperature in the Arctic is 5 degrees
C for a year, the simple model predicts 77 Pg of carbon will be released 
in the form of CO2 whilst the integral model says 35 Pg will be released. 
Neither of these values are particularly accurate, but through the use of 
the proposed model we have significantly less emissions and we are able to 
capture the effect that temperatures decrease as we move pole-ward and 
vice-versa as we move towards the equator. 

Next steps include modifying the insolation coefficient distribution presented
in McGehee to better reflect how our temperature scaling factors should 
behave as a function of latitude. 
%}
%% Visualizing coefficient distribution 
coeff = @(x) 1 + ((5*pi/12)-x)*((10/273.15)/(pi/12));

x = linspace(pi/3, pi/2, 50);
for i = 1:50
    y(i) = coeff(x(i));
end 

plot(x,y)
xlim([pi/3 pi/2])
ylim([0.8 1.1])
yline(1, '-.r')
xlabel('Latitude')
ylabel('Scaling factor')
%% FUNCTIONS 
function m = mthaw(tau, m)
%{
Given a total mass of carbon within a surface slab and the temperature at
the surface we compute the estimated mass of carbon at above freezing
temperatures. 

Args:  tau - non-dimensional temperature
       m - mass of carbon in surface slab

Returns: m - mass of carbon susceptible to thaw

I use the equilibrium temperature gradient function which is 
defined by the boundary conditions T_{z=0} = tau_s, T_{z=Zb} = Tb, 
(our original temperature gradient).  This will change in the future, 
but it was easier to capture the idea of a temperature gradient in this way
than through the use of our current idea.  We just need to illustrate 
the concept that susceptible thaw mass changes with surface temperature.
Also our new idea would require imbedding the mthaw function within the
EBM which detracts from the main point of this script. 
%}

T = tau.*273.15-273.15; %non-dimensional temperature to degrees C
zd = 3; %surface slab thickness
rho = m./zd; %carbon per m^3 within surface slab (slab is unitary in area)
Td = -2.35; %constant temperature at bottom of slab in degrees C

if T>=0
    z0 = (T.*zd)./(T-Td);   %solution for zero temp depth 
else 
    z0 = 0; %depth is zero if surface temp is below freezing 
end 

m = rho.*z0;  %thaw mass in an area unit of Arctic soil 
end 
function r2 = Q10(r1, t1, t2, sensitivity)
%{
Q10 sensitivity function 
Args:  r1 - reference decomposition rate from literature
       t1 - temperature associated with r1
       t2 - new temperature 
       sensitivity - Q10 value 
%}
    r2 = r1*sensitivity.^((t2-t1)./10);
end
function r = rateCO2(tau)
%{
Function for scaling rate of decomposition according to Q10 sensitivity
Args: tau - nondimensional temperature 
Returns: r - fraction of carbon released within a year as CO2

uses pool speeds, soil type and aerobic conditions 
%}

%Constants
gammaA = 0.025; %fraction of carbon in active pool
gammaS = 0.45; %fraction of carbon in slow pool
fms = 0.57; %mass fraction of carbon in mineral soils
fo = 0.22; %mass fraction of carbon in organic soils 
fa = 0.07; %mass fraction of carbon in arctic river deltas
fy = 0.14; %mass fraction of carbon in yedoma
sigma_yms = 0.79; %fraction of yedoma deposits in mineral soils 
sigma_yo = 0.21; %fraction of yedoma deposits in organic soils
omega_yan = 0.33; %fraction of anaerobic decomposition in yedoma 
A_msan = 0.05; %Area fraction in mineral with anaerobic decomp
A_oan = 0.8; %Area fraction in organic with anaerobic decomp
chi_ms = 0.25; %oxidation fraction mineral under anaerobic
chi_o = 0.6; %oxidation fraction organic under anaerobic
R_ana = 0.1; %ratio of decomp speeds anaerobic:aerobic 
C_tot = 1300*10.^12; %total carbon stocks in permafrost 
Q10_a = 1.5; %Q10 sensitivity aerobic
Q10_an = 3; %Q10 sensitivity anaerobic 
k_a = 1./0.37; %decomp rate active at 5 degrees C (1/MRT)
k_s = 1./7.46; %decomp rate slow  at 5 degrees C
k_e = 0.01; %decomp ebullition ad hoc for now 

%define rate scaling functions
T = tau.*273.15-273.15;
ratescale_a = @(T) Q10(k_a, 5, T, Q10_a);
ratescale_s = @(T) Q10(k_s, 5, T, Q10_a);

%logic from "C Release" document 
rate_co2 = (gammaA*ratescale_a(T) + gammaS.*ratescale_s(T)).*((fms.*(1-A_msan)+ fo.*(1-A_oan) + fy.*(1-omega_yan)) + R_ana.*(fms.*A_msan.*chi_ms + fo.*A_oan.*chi_o));

r = rate_co2;
end
function r = rateCH4(tau)
%{
Function for scaling rate of decomposition according to Q10 sensitivity
Args: tau - nondimensional temperature 
Returns: r - fraction of carbon released within a year as CH4

uses pool speeds, soil type and aerobic conditions 
%}

%Constants
gammaA = 0.025; %fraction of carbon in active pool
gammaS = 0.45; %fraction of carbon in slow pool
fms = 0.57; %mass fraction of carbon in mineral soils
fo = 0.22; %mass fraction of carbon in organic soils 
fa = 0.07; %mass fraction of carbon in arctic river deltas
fy = 0.14; %mass fraction of carbon in yedoma
sigma_yms = 0.79; %fraction of yedoma deposits in mineral soils 
sigma_yo = 0.21; %fraction of yedoma deposits in organic soils
omega_yan = 0.33; %fraction of anaerobic decomposition in yedoma 
A_msan = 0.05; %Area fraction in mineral with anaerobic decomp
A_oan = 0.8; %Area fraction in organic with anaerobic decomp
chi_ms = 0.25; %oxidation fraction mineral under anaerobic
chi_o = 0.6; %oxidation fraction organic under anaerobic
R_ana = 0.1; %ratio of decomp speeds anaerobic:aerobic 
C_tot = 1300*10.^12; %total carbon stocks in permafrost 
Q10_a = 1.5; %Q10 sensitivity aerobic
Q10_an = 3; %Q10 sensitivity anaerobic 
k_a = 1./0.37; %decomp rate active at 5 degrees C (1/MRT)
k_s = 1./7.46; %decomp rate slow  at 5 degrees C
k_e = 0.01; %decomp ebullition ad hoc for now 

%define rate scaling functions
T = tau.*273.15-273.15;
ratescale_a = @(T) Q10(k_a, 5, T, Q10_a);
ratescale_s = @(T) Q10(k_s, 5, T, Q10_a);

%logic from "C Release" document 
rate_ch4 = R_ana.*(gammaA*ratescale_a(T) + gammaS.*ratescale_s(T)).*(fms.*A_msan.*(1-chi_ms)+ fo.*A_oan.*(1-chi_o))+ k_e.*(fa + fy.*omega_yan);

r = rate_ch4;
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
function m = integralmethod(tau, mass)
%{
implements temperature distribution mechanism for determining thaw mass 
Args: tau - non-dimensional temperature at surface
      mass - mass of carbon in a unit area 
Returns: m - tuple of total mass of carbon respired as [mCO2 mCH4]

temperature coefficient function was decided to be linear for
the time being, with a central value of 1 located at the midpoint
between pi/3 and pi/2.  Slope is such that total temperature difference between
60N and the pole is 15 degrees C. 
%}

%Linear function for temperature scaling factor 
coeff = @(x) 1 + ((5*pi/12)-x)*((7.5/273.15)/(pi/12));

mrCO2 = @(y,t,mass) mthaw(t*coeff(y), mass).*rateCO2(t*coeff(y)); %mass of carbon respired as CO2 in a year
mrCH4 = @(y,t,mass) mthaw(t*coeff(y), mass).*rateCH4(t*coeff(y)); %mass of carbon respired as CH4 in a year 

%integrate over Arctic 
m(1) = integral(@(y) mrCO2(y,tau,mass).*2.*pi.*(6371.*1000).^2.*cos(y).*cos(y), pi./3, pi./2, 'ArrayValued', true);
m(2) = integral(@(y) mrCH4(y,tau,mass).*2.*pi.*(6371.*1000).^2.*cos(y).*cos(y), pi./3, pi./2, 'ArrayValued', true);
end