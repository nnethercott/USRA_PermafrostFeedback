%% INTEGRAL METHOD 
%initial params 
m = 1095.*10.^12;
A = integral(@(x) (6371.*1000).^2.*2.*pi.*sin(x), 0, pi/6);
Z_L = 15;
cdensity = @(m) m./(A.*Z_L); 
options = optimset('Display','off');
x0 = 0.9; %solver x0

F_A = 104;
Q = 201.73; 
F_O = 9.75;
delta = 0.6;

%arrays 
ch4conc(1) = CH485data(1);    
co2conc(1) = CO285data(1);
co2_loss(1) = 0;
ch4_loss(1) = 0;

%FUNCTIONS
%define latitude dependent insolation function from McGehee
relinsol = @(x,y) 2.*sqrt(1 - (sqrt(1 - sin(y).^2).*sin(0.4101524).*cos(x) - sin(y).*cos(0.4101524)).^2)/pi.^2;
coeff = @(y) integral(@(x) relinsol(x,y), 0, 2.*pi);
insol = @(y) 343.*coeff(y);

for i = 1:3
    m = m - co2_loss(i) - ch4_loss(i);
    rhoC = cdensity(m);
    
    %equilibrium states
    eqsf = @(y) fsolve(@(tau) eqtemp(tau, F_O, F_A, insol(y), co2conc(i), delta, ch4conc(i)),x0, options);
    eqaf = @(y) eqairtemp(eqsf(y), F_O, F_A, insol(y), co2conc(i), delta, ch4conc(i));
    Z = @(y) fsolve(@(z) eqtempgrad(eqsf(y), F_O, F_A, insol(y), co2conc(i), delta, ch4conc(i), z)-1,0, options);

    %functions for thaw mass at a latitude
    mCO2 = @(y) integral(@(z) rateCO2(eqtempgrad(eqsf(y), F_O, F_A, insol(y), co2conc(i), delta, ch4conc(i), z)).*rhoC, Z(y), 0,'ArrayValued', true);
    mCH4 = @(y) integral(@(z) rateCH4(eqtempgrad(eqsf(y), F_O, F_A, insol(y), co2conc(i), delta, ch4conc(i), z)).*rhoC, Z(y), 0,'ArrayValued', true);

    %compute 
    m_totalCO2_2 = integral(@(y) 2.*pi.*(6371.*1000).^2.*cos(y).*mCO2(y), pi/3, pi/2);
    m_totalCH4_2 = integral(@(y) 2.*pi.*(6371.*1000).^2.*cos(y).*mCH4(y), pi/3, pi/2);
    
    %my version of uptake integral
    rads=linspace(pi/3,pi/2,300);
    dy=(pi/2-pi/3)/300;
    sink2=0;
    for j=1:300
        sink2=sink2+2.*pi.*(6371.*1000).^2.*cos(rads(j)).*npp(co2conc(i),eqaf(rads(j)))*dy;
    end 

    m_netCO2_2 = m_totalCO2_2-sink2;
    m_netCH4_2 = m_totalCH4_2;
    
    co2_loss(i+1) = m_netCO2_2;
    ch4_loss(i+1) = m_netCH4_2;

    %diverge from rcp concentrations
    new_concentrations = concentration([sum(co2_loss) sum(ch4_loss)], [CO285data(i+1) CH485data(i+1)]);
    co2conc(i+1) = new_concentrations(1);
    ch4conc(i+1) = new_concentrations(2);
    
    eqnord = eqsf(pi/2)
    eqsud = eqsf(pi/3)
    eqanord = eqaf(pi/2)
    eqasud = eqaf(pi/3)
    znord = Z(pi/2)
    zsud = Z(pi/3)
    sink2
    m_netCO2_2
    m_netCH4_2
    m
end 
