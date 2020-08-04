The goal of this repo is to implement permafrost dynamics within an existing energy balance model to capture non-linear feedback terms
which occur as a consequence of increasing Arctic temperatures. 

All work is currently on the dev branch.

Work presented in this repo is built off [An energy balance model for paleoclimate transitions](https://doi.org/10.5194/cp-15-493-2019) and is a joint effort with the Department of Mathematics and Statistics at the University of Guelph.

General
=======

... The surface slab is split into two components; ocean and land. The
temperature of the ocean section is assumed to be uniform while the
temperature of the land portion changes with depth. At the surface of
both mediums their temperatures are equal, i.e. $T_{O}\equiv T_{S}(0,t)$
where $T_{O}$ is the temperature of the ocean and $T_{S}(z,t)$
represents the temperature of the soil at time $t$ and at depth
$-Z_{L} \leq z \leq 0$. Land and ocean areas are assumed to be
“well-mixed”. We assume fluxes within the soil section to be
one-dimensional, acting normal to the surface-atmosphere interface.
Applying energy conservation in the ocean layer and the modelling the
heat equation in the soil yields

$$\rho_{O}c_{O}Z_{O}\frac{\partial T_{S}(0,t)}{\partial t} = -K_{L}\frac{\partial T_{S}(0,t)}{\partial z} + q(0,t)$$

$$\rho_{L}c_{L}\frac{\partial T_{S}(z,t)}{\partial t} = K_{L}\frac{\partial^{2}T_{S}(z,t)}{\partial z^{2}}$$

$$T_{S}(-Z_{L},t) = T_{Z_{L}}$$

where $\rho$ is density in kg m$^{-3}$, $c$ is specific heat capacity in
J kg$^{-1}$ K$^{-1}$, $Z_{O}$ is the depth of our ocean slab in m,
$K_{L}$ the thermal conductivity of land in W m$^{-1}$ K$^{-1}$, and $q$
is the total flux source at the surface in W m$^{-2}$. Note that by
equation (3) we impose the boundary condition that there exists a depth
where temperature remains constant in the soil. Substituting our EBM
parameters, the system of equations governing the climate model become

$$\begin{aligned}
\rho_{A}c_{A}Z_{A}\frac{dT_{A}(t)}{dt} & =  F_{A} + F_{C} +\xi_{A}Q + \eta I_{S} - I_{A}\\
 \rho_{O}c_{O}Z_{O}\frac{\partial T_{S}(0,t)}{\partial t} & =  -K_{L}\frac{\partial T_{S}(0,t)}{\partial z} + F_{O} - F_{C} + (1-\alpha)F_{S} - I_{S} + \beta I_{A}\\
 \rho_{L}c_{L}\frac{\partial T_{S}(z,t)}{\partial t} & = K_{L}\frac{\partial^{2}T_{S}(z,t)}{\partial z^{2}}\\
 T_{S}(-Z_{L},t) & = T_{Z_{L}}\end{aligned}$$

We non-dimensionalize the above system using the following dimensionless
variables

$$\begin{aligned}
s = \frac{\sigma T_{R}^{4}}{c_{O}\rho_{O}Z_{O}T_{R}}t, \qquad
\tau_{S} = \frac{T_{S}}{T_{R}}, \qquad \tau_{A} = \frac{T_{A}}{T_{R}}, \qquad \zeta = \frac{z}{Z_{L}},\end{aligned}$$

where $\sigma$ is the Stefan-Boltzmann constant and $T_{R}=273.15 K$ is
our reference temperature. Here one unit of dimensionless time $s$
corresponds to the amount of actual time required for the earth at
$T_{R}$ to radiate an equivalent amount of energy which is stored in our
ocean slab at $T_{R}$. The non-dimensional system becomes

$$H_{1}\frac{d\tau_{A}(s)}{ds} =  f_{A} + f_{C} +\xi_{A}q + \eta i_{S} - i_{A}$$

$$\frac{\partial \tau_{S}(0,s)}{\partial s} + H_{2}\frac{\partial \tau_{S}(0,s)}{\partial \zeta} = f_{O} - f_{C} + (1-\alpha)f_{S} - i_{S} + \beta i_{A}$$

$$H_{3}\frac{\partial \tau_{S}(\zeta,s)}{\partial s} = \frac{\partial^{2}\tau_{S}(\zeta,s)}{\partial \zeta^{2}}$$

$$\tau_{S}(-1, s) = \tau_{Z_{L}}$$

where

$$\begin{aligned}
    H_{1} = \frac{\rho_{A}c_{A}Z_{A}}{\rho_{O}c_{O}Z_{O}}, & \qquad H_{2} = \frac{K_{L}}{Z_{L} \sigma T_{R}^{3}}, \qquad H_{3} = \frac{\rho_{L}c_{L}Z_{L}^{2}\sigma T_{R}^{3}}{K_{L}\rho_{O}c_{O}Z_{O}}, \qquad \tau_{Z_{L}} = \frac{T_{Z_{L}}}{T_{R}}, \\
    \qquad f_{A} = \frac{F_{A}}{\sigma T_{R}^4}, & \qquad f_{C} = \frac{F_{C}}{\sigma T_{R}^4}, \qquad f_{O} = \frac{F_{O}}{\sigma T_{R}^4}, \qquad q = \frac{Q}{\sigma T_{R}^4}, \qquad i_{S} = \frac{I_{S}}{\sigma T_{R}^4}, \\
    i_{A} = \frac{I_{A}}{\sigma T_{R}^4} &\\\end{aligned}$$

We are interested in equilibrium solutions to our system of equations
where $\tau_{S}(\zeta,s) = \tau_{S}(\zeta)$ and
$\tau_{A}(s) = \tau_{A}$. Our system thus becomes

$$0 =  f_{A} + f_{C} +\xi_{A}q + \eta i_{S} - i_{A}$$

$$H_{2}\frac{\partial \tau_{S}(0)}{\partial \zeta} = f_{O} - f_{C} + (1-\alpha)f_{S} - i_{S} + \beta i_{A}$$

$$0 = \frac{\partial^{2}\tau_{S}(\zeta)}{\partial \zeta^{2}}$$

$$\tau_{S}(-1) = \tau_{Z_{L}}$$

Solving this system we obtain

$$\tau_{S}(\zeta) = \tau_{Z_{L}} + \frac{\zeta + 1}{H_{2}}[f_{O} -(1-\beta)f_{C} + (1-\alpha)f_{S} -(1-\beta\eta)i_{S} + \beta( f_{A} + \xi_{A}q)]$$

Permafrost
==========

Background
----------

Perennially frozen soils contain massive amounts of carbon which become
increasingly susceptible to thaw as global temperatures rise. Carbon
released from the permafrost contributes to atmospheric GHG burdens
which in turn accelerates the warming process and results in further
carbon losses from these soils. This is permafrost feedback. In this
permafrost module we limit our analysis to carbon stored within the top
3 m of High-Arctic permafrost soils, where the Arctic is defined as any
region with latitude greater than or equal to 60N. Carbon emission rates
are determined from aerobic conditions, soil type and pool decomposition
speed. We consider two soil types described by their organic carbon
content; mineral soils which are $<$20%C by weight and organic soils
which are $\geq$20%C by weight (@athaw). We use a 3-pool decomposition
model where carbon falls into a passive, slow or active decomposition
pool based on its lability. Carbon in the active pool is respired by
permafrost soils within a year, carbon in the slow pool on a decadal
scale, and carbon in the passive pool over hundreds to thousands of
years (@schadel). We assume carbon in the passive pool is inert. In the
simplest case, carbon is assumed to be respired from the permafrost in
the form of CO$_{2}$ or CH$_{4}$. What determines this is the
availability of oxygen in the decay process. It is assumed that aerobic
decomposition produces strictly CO$_{2}$ while anaerobic decomposition
produces CH$_{4}$. Of the carbon which decomposes anaerobically, a
fraction is oxidized before reaching the surface-atmosphere interface.
We assume this oxidized fraction is released as CO$_{2}$. Both CO$_{2}$
and CH$_{4}$ are “well-mixed” within the atmosphere. The atmospheric
lifetime of CH$_{4}$ is roughly a decade due to removal by the OH
radical. We assume CH$_{4}$ is broken down into CO$_{2}$ during this
process and that after 10 years the entirety of a CH$_{4}$ sample will
have been converted to CO$_{2}$.

Decomposition
-------------

Suppose carbon is distributed within the soil according to some function
$m(z,t)$. Using a temperature-dependent decay rate constant, $k(T)$, we
express the decomposition of carbon at a specific depth in the soil as

$$\frac{dm(z,t)}{dt} = -k(T_{S}(z,t))m(z,t)$$

Let us denote by $m_{O}(z,t)$ the mass distribution of carbon in pools
where aerobic decomposition takes place, $m_{X}(z,t)$ the mass
distribution of carbon in anaerobic pools where oxidation occurs, and
$m_{N}(z,t)$ the mass distribution of carbon in anaerobic pools with no
oxidation. These pools decompose according to equation (13) with their
own decay rate constants. We then express the decay of permafrost carbon
into CO$_{2}$ and CH$_{4}$ as follows

$$\begin{split}
    \frac{dCO_{2}(t)}{dt} = \frac{M_{CO_{2}}}{M_{C}}\int_{Z_{L}}^{0} k_{O}(T_{S}(z,t))m_{O}(z,t) + k_{X}(T_{S}(z,t))m_{X}(z,t) dz
    \end{split}$$

$$\frac{dCH_{4}(t)}{dt} = \frac{M_{CH_{4}}}{M_{C}}\int_{Z_{L}}^{0} k_{N}(T_{S}(z,t))m_{N}(z,t)dz$$

Here $M_{C}$, $M_{CO_{2}}$ and $M_{CH_{4}}$ are the molar masses of C,
CO$_{2}$ and CH$_{4}$. Equations (14) and (15) describe release rates of
each gas. One should note that since we’ve assumed carbon is only in the
top 3 m, the lower bound on the integrals above can be replaced with
$-3$. Once in the atmosphere we must consider sink terms introduced by
the OH radical which acts to remove CH$_{4}$ and convert it to CO$_{2}$
(by assumption). Modifying the above equations we obtain

$$\begin{split}
    \frac{dCO_{2}(t)}{dt} = \frac{M_{CO_{2}}}{M_{C}}\int_{Z_{L}}^{0} k_{O}(T_{S}(z,t))m_{O}(z,t) + k_{X}(T_{S}(z,t))m_{X}(z,t) dz
    \\ + (\frac{M_{CO_{2}}}{M_{CH_{4}}})r_{CH_{4}}CH_{4}(t)
    \end{split}$$

$$\frac{dCH_{4}(t)}{dt} = \frac{M_{CH_{4}}}{M_{C}}\int_{Z_{L}}^{0} k_{N}(T_{S}(z,t))m_{N}(z,t)dz - r_{CH_{4}}CH_{4}(t)$$

where $r_{CH_{4}}$ is the decay rate constant for the removal of
CH$_{4}$ by OH.

We now wish to model how the atmospheric concentrations of
permafrost-borne CO$_{2}$ and CH$_{4}$ change with time while
considering sink terms. Let the total mass of the atmosphere be given by
$m_{A}(t)$, we have that

$$\frac{d\mu(t)}{dt} = (\frac{M_{A}}{M_{CO_{2}}})\frac{\dot{CO_{2}(t)}}{m_{A}(t)+\dot{CO_{2}(t)}+\dot{CH_{4}(t)}}\times 10^{6} + (\frac{M_{CO_{2}}}{M_{CH_{4}}})r_{CH_{4}}\nu(t)$$

$$\frac{d\nu(t)}{dt} = (\frac{M_{A}}{M_{CH_{4}}}) \frac{\dot{CH_{4}(t)}}{m_{A}(t)+\dot{CO_{2}(t)}+\dot{CH_{4}(t)}}\times 10^{6} -r_{CH_{4}}\nu(t)$$

where $M_{A}$ denotes the molar mass of the atmosphere and $r_{CH_{4}}$
the atmospheric decay rate of CH$_{4}$ due to removal by OH. Here
$\dot{\mu}$ and $\dot{\nu}$ are in ppm $t^{-1}$. Note the emissions of
CO$_{2}$ and CH$_{4}$ increase the total mass of the atmosphere in the
instant they are added so we account for this by including their effects
in the denominator of the mass ratios. Solutions to equations (16) and
(17) allow us to determine the additional permafrost contributions to
the global atmospheric GHG burden in a given year.

Modelling
=========

We first construct a discrete version our mass distribution function by
partitioning the top 3 m of our surface slab into $P$ pools and dividing
the total permafrost carbon amongst them. Denote by
$m^{i}(t)=m(z^{i}, t)$ the mass distribution of carbon in pool $i$ with
$z^{i} = -3(\frac{i-1}{P})$. Likewise, define by $T^{i}(t)=T(z^{i},t)$
the temperature of the i$^{th}$ pool. In order to determine the effect
of permafrost feedback on the EBM, we need solutions to equations (16)
and (17) which in turn requires knowledge of $m(z,t)$. In the discrete
case, the system of equations governing this problem becomes

$$\begin{split}
    \frac{dm^{i}(t)}{dt} &= -k_{O}(T_{S}^{i}(t))m^{i}_{O}(t)-k_{X}(T_{S}^{i}(t))m^{i}_{X}(t)-k_{N}(T_{S}^{i}(t))m^{i}_{N}(t)\\
    \frac{dCO_{2}(t)}{dt} &= \frac{M_{CO_{2}}}{M_{C}}\sum_{i=1}^{P} k_{O}(T_{S}^{i}(t))m^{i}_{O}(t) + k_{X}(T_{S}^{i}(t))m^{i}_{X}(t) +  (\frac{M_{CO_{2}}}{M_{CH_{4}}})r_{CH_{4}}CH_{4}(t)\\
    \frac{dCH_{4}(t)}{dt} &= \frac{M_{CH_{4}}}{M_{C}}\sum_{i=1}^{P} k_{N}(T_{S}^{i}(t))m^{i}_{N}(t) - r_{CH_{4}}CH_{4}(t)\\
    \end{split}$$

The solution to this system of $P$+2 differential equations is a
function $u(t)$ such that

$$\begin{aligned}
    \dot{u(t)} &= 
    \begin{bmatrix}
           \frac{m^{1}(t)}{dt} \\
           \frac{m^{2}(t)}{dt} \\
           \vdots \\
           \frac{m^{P}(t)}{dt} \\
           \frac{dCO_{2}(t)}{dt} \\
           \frac{dCH_{4}(t)}{dt} \\
         \end{bmatrix},\:
         u(t_{0}) = \begin{bmatrix}
           m^{1}(t$_{0}$) \\
           m^{2}(t$_{0}$) \\
           \vdots \\
           m^{P}(t$_{0}$) \\
           0 \\
           0 \\
           \end{bmatrix}
  \end{aligned}$$

Here we assume that no permafrost carbon is present in the atmosphere in
the year $t_{0}$. We numerically solve for $u(t)$ using a
Adams-Bashforth 3-step method along different RCP projections with
$u(t_{1})$ and $u(t_{2})$ being obtained from a third order Runge-Kutta
method.

Results
=======

When considering the impact of including permafrost terms in a climate
model, the natural assumption is that the onset of a bifurcation (if
any) will be expedited. In the very least we assume additional carbon
burdens introduced by permafrost thaw will result in proportionally
warmer Arctic surface temperatures. It would be nice if numerical
analysis validated our intuition. Using RCP projections and their
extensions (ECPs) in tandem with historical records we are able to
prescribe realistic atmospheric carbon and methane concentrations on an
annual basis from 1765 until 2500. We will limit our simulations to the
range of years between 1900 and 2300 inclusive. When plotting the
permafrost and no-permafrost situations we see our initial inklings are
supported. In Fig.1 we find that the inclusion of permafrost terms
results in a bifurcation along RCP 8.5 in the year 2107 whereas this
occurs without permafrost considerations over a decade later in 2118.
One should also note the divergence of RCP 6.0 simulated with permafrost
feedback from its reference curve.

![Arctic surface temperatures simulated along the four RCP pathways from
1900 to 2300 while considering permafrost feedback.](RCPs.jpg "fig:")
[fig:my~l~abel]

By the year 2299 the entirety of the initial carbon stocks in permafrost
soils - 1095 PgC - have been lost along RCP 8.5. This estimate falls
well outside the range outlined in @mcguire which reports values of 641
Pg C losses to 167 Pg C gains on the extremes. This most likely is a
consequence of using steady state methods to model the heat equation in
the surface slab as well as due to simplifications made in our modelling
efforts (uniform surface slab, no latent heat, 1-D heat equation, etc.).
We also do not model NPP and other plant effects which might function to
offset these HR losses. For RCPs 4.5 and 2.0 there is no observed thaw.
In Fig. 1 we’ve prescribed the EBM with intermediate Arctic atmospheric
and ocean heat fluxes of $F_{A} = 104$ W m$^{-2}$ and $F_{O} = 10$ W
m$^{-2}$ respectively. By instead running our simulation using a linear
interpolation of flux values as reported by @koenigk (Table 1) we come
across an interesting result. For sufficiently large annual temperature
fluctuations, $\Delta T_{ann}$, the inclusion of permafrost feedback
along RCP 6.0 under these conditions causes a bifurcation to the warm
solution branch to occur where such a feature does not otherwise exist.
For appropriately small values of $\Delta T_{ann}$, divergence from the
baseline curve is still observed.

[H]

[H]<span>m<span>1.2cm</span>m<span>2.0cm</span>m<span>2.2cm</span>m<span>2.2cm</span></span>\
**<span>Year</span> & **<span>Scenario</span> & **<span>$F_{A}$ (W
m$^{-2}$)</span> & **<span>$F_{O}$ (W m$^{-2}$)</span>\
1850 & Historical & 107.28 & 9.75\
2100 & RCP 2.6 & 104.03 & 13.00\
2100 & RCP 8.5 & 97.53 & 19.5\
********

[tab:my~l~abel]

<span>cc</span> ![image](kbRCPs12.jpg) & ![image](kbRCPs11.jpg)\
(a) T$_{ann}$ = 12C& (b) T$_{ann}$ = 11C\
\

When simulating the EBM using instead $F_{A}$ projections obtained from
@yang we find that a bifurcation results along RCP 6.0 regardless of the
inclusion of permafrost feedback. In this scenario $F_{A}$ increases
with time (unlike in the previous case) and reaches a peak value of 129
W m$^{-2}$ by 2100. Note also that a bifurcation takes place along RCP
4.5 when one considers permafrost feedback.

![Arctic surface temperatures simulated along RCPs from 1900 to 2300
using @yang $F_{A}$ data. $F_{O}$ is as before, and both fluxes increase
linearly until 2100 after which they are held
constant.](yangRCPs.jpg "fig:") [fig:my~l~abel]

By changing the value of Q$_{10}$ in our temperature-dependent decay
rate equations (see section 5.2) we are modifying the sensitivity of our
decomposition processes to soil temperatures. Of course this change is
only observed where thaw occurs. In Fig. 4 we show plots for RCPs 8.5
and 6.0 with varying Q$_{10}$ values. We revert back to using the same
intermediate values of 10 W m$^{-2}$ and 104 W m$^{-2}$ for $F_{O}$ and
$F_{A}$ respectively in these scenarios.

<span>cc</span>\
\

It was seen earlier that RCP 8.5 exhibits a bifurcation regardless of
permafrost considerations, with the inclusion of permafrost feedback
speeding up the process by a decade. In Fig. 4a we see that modifying
Q$_{10}$ has very little impact on the timing of our bifurcation. From
these two pieces of information we are left with the conclusion that the
bifurcation along RCP 8.5 is largely a consequence of other forcing
variables and that the inclusion of permafrost terms does not result in
any characteristic change. Under intermediate conditions, changing the
Q$_{10}$ parameter is not sufficient on its own to cause a bifurcation
along RCP 6.0. In Fig. 4b we see that increasing the Q$_{10}$ actually
results in decreasing surface temperatures – a relationship which goes
against our intuition but is easily explained. Q$_{10}$ sensitivity is
described with respect to a reference rate at a reference temperature
which in our case were decomposition rates obtained from @schadel for
aerobic incubations carried out at 5C. Since soil temperatures are below
freezing in the situation above, increasing the sensitivity of our decay
rates to temperature has the effect of reducing decomposition rates at
cooler temperatures.

Supplementary
=============

Physical parameters
-------------------

[H]<span>m<span>1.2cm</span>m<span>5.9cm</span>m<span>2.2cm</span>m<span>2cm</span></span>
**Symbol & **Description & **Value & **Reference\
$C_{tot}$ & Total C stocks & 1300 PgC & @hugelius\
$S_{A}$ & Surface area of Arctic region & 3.4210$^{13}$ m$^{2}$ & N/A\
$S_{P}$ & Permafrost surface area & 1.7810$^{13}$ m$^{2}$ & @hugelius\
$Z_{L}$ & Soil slab depth & 15 m & @gtnp\
$T_{S}(Z_{L},t)$ & Constant temperature at bottom of surface slab &
-3.054C & @gtnp\
$K_{L}$ & Thermal conductivity of land & 2 W m$^{-1}$ K$^{-1}$ &
**<span>SOURCE</span>\
$\gamma_{a,ms}$ & Fraction of carbon available for decomposition in
active pool, mineral soils & 0.013 & @schadel\
$\gamma_{s,ms}$ & Fraction of carbon available for decomposition in slow
pool, mineral soils & 0.107 & @schadel\
$\gamma_{a,o}$ & Fraction of carbon available for decomposition in
active pool, organic soils & 0.015 & @schadel\
$\gamma_{s,o}$ & Fraction of carbon available for decomposition in slow
pool, organic soils & 0.293 & @schadel\
*f$_{ms}$* & Mass fraction of carbon in mineral soils & 0.708 &
@hugelius\
*f$_{o}$* & Mass fraction of carbon in organic soils & 0.292 &
@hugelius\
*A$_{ms,an}$* & Area fraction of mineral soil with anaerobic
decomposition & 0.05 & @schneider\
*A$_{o,an}$* & Area fraction of organic soils with anaerobic
decomposition & 0.8 & @frolking\
$_{ms}$ & Fraction of CH$_{4}$ oxidized in mineral soils & 0.25 &
@schneider\
$_{o}$ & Fraction of CH$_{4}$ oxidized in organic soils & 0.6 &
@schneider\
R$_{an/a}$ & Ratio of decomposition rates in aerobic and anaerobic
conditions & 0.1 & @schneider\
$Q_{10}$ & Q$_{10}$ temperature sensitivity & 2.5 & @schadel\
$\kappa_{a,ms}$ & Active pool decomposition rate constant at 5 C in
mineral soils under aerobic conditions& 2.8986 & @schadel\
$\kappa_{s,ms}$ & Slow pool decomposition rate constant at 5 C in
mineral soils under aerobic conditions& 0.1318 & @schadel\
$\kappa_{a,o}$ & Active pool decomposition rate constant at 5 C in
organic soils under aerobic conditions& 2.4390 & @schadel\
$\kappa_{s,o}$ & Slow pool decomposition rate constant at 5 C in organic
soils under aerobic conditions& 0.1387 & @schadel\
$r_{CH_{4}}$ & Atmospheric decay rate of CH$_{4}$ & 0.0762 &
**<span>SOURCE</span>\
& & &\
 ************

Derivation of $\dot{m_{O}}$, $\dot{m_{X}}$, and $\dot{m_{N}}$
-------------------------------------------------------------

Suppose the carbon mass distribution in permafrost soils is given by
$m(z,t)$. We assume carbon pools corresponding to anaerobic, anaerobic
with oxidation and anaerobic without oxidation are well-mixed such that
the mass of carbon in each pool can be expressed in the form fraction
times $m(z,t)$. Our goal is to construct these fractions. Using
assumptions outlined in section 2.1 as well as parameters listed in
Table 1 we arrive at

$$m_{O}(z,t) = m(z,t)(f_{ms}(1-A_{ms,an}) + f_{o}(1-A_{o,an}))$$

$$m_{X}(z,t) = m(z,t)(f_{ms}A_{ms,an}\chi_{ms} + f_{o}A_{o,an}\chi_{o})$$

$$m_{N}(z,t) = m(z,t)(f_{ms}A_{ms,an}(1-\chi_{ms}) + f_{o}A_{o,an}(1-\chi_{o}))$$

Now by equation (4) all we need to convert the above equations into
emission rates is the introduction of a decomposition rate constant,
$k$. We express this value as a sum of active pool and slow pool
contributions in mineral and organic soils.

$$\begin{split}
    \frac{dm_{O}(z,t)}{dt} = -m(z,t)(f_{ms}(1-A_{ms,an})[\gamma_{a,ms}k_{a,ms}(T_{S}(z,t)) \\+ \gamma_{s,ms}k_{s,ms}(T_{S}(z,t))] + f_{o}(1-A_{o,an})[\gamma_{a,o}k_{a,o}(T_{S}(z,t)) \\
    + \gamma_{s,o}k_{s,o}(T_{S}(z,t))])
    \end{split}$$

$$\begin{split}
    \frac{dm_{X}(z,t)}{dt} = -m(z,t)(f_{ms}A_{ms,an}\chi_{ms}[\gamma_{a,ms}R_{an/a}k_{a,ms}(T_{S}(z,t)) \\+ \gamma_{s,ms}R_{an/a}k_{s,ms}(T_{S}(z,t))] + f_{o}A_{o,an}\chi_{o}[\gamma_{a,o}R_{an/a}k_{a,o}(T_{S}(z,t)) \\+ \gamma_{s,o}R_{an/a}k_{s,o}(T_{S}(z,t))])
    \end{split}$$

$$\begin{split}
    \frac{dm_{N}(z,t)}{dt} = -m(z,t)(f_{ms}A_{ms,an}(1-\chi_{ms})[\gamma_{a,ms}R_{an/a}k_{a,ms}(T_{S}(z,t)) \\+ \gamma_{s,ms}R_{an/a}k_{s,ms}(T_{S}(z,t))] + f_{o}A_{o,an}(1-\chi_{o})[\gamma_{a,o}R_{an/a}k_{a,o}(T_{S}(z,t)) \\+ \gamma_{s,o}R_{an/a}k_{s,o}(T_{S}(z,t))])
    \end{split}$$

here $k(T_{S}(z,t))$ is an annual effective decomposition rate given by

$$k(T_{S}(z,t)) = \frac{1}{2\pi}\int_{0}^{2\pi}g(T,t)dt$$

where

$$g(T,t) = 
    \begin{cases}
    \kappa\cdot Q_{10}^{\frac{T_{S}(z,t)-5}{10}} & T + \Delta T \sin{t}>0 \\
    0   & else
    \end{cases}$$

With equations (16) and (17) we capture the effects of seasonal
temperature fluctuations so that even when the annual mean temperature
is below freezing there can still be permafrost thaw as a consequence of
higher temperatures experienced in the spring and summer months. By
substituting equations (12) through (14) into equations (5) and (6) we
are able to describe the emission rates of CO$_{2}$ and CH$_{4}$.

Absorptivity CH$_{4}$
---------------------

In @reto we see that the authors attribute 50% of the modern day
absorption factor to water vapour, 25% to clouds, 19% to CO$_{2}$ and 7%
to other gasses. Of the other gasses, CH$_{4}$ accounts for 1 seventh of
the their absorption, so CH$_{4}$ contributes 1% overall. We will lump
the remaining effects into the fraction for CO$_{2}$. The the modern-day
value for $\eta$ is 0.9, so

$$0.9 = 1-(1-50x)(1-25x)(1-24x)(1-x)$$

Solving the quartic for $x$ yields $x = 0.01494$. We also have that

$$x = 1 - e^{-\nu\cdot \frac{M_{CH_{4}}}{M}\times10^{-9}k_{CH_{4}}\frac{P_{A}}{g}}$$

Taking the molar mass of methane to be 16.04 g mol$^{-1}$, the molar
mass of atmosphere to be 28.97 g mol$^{-1}$ (so mmCH$_{4}$/mmAtm =
0.554) and year 2000 levels of CH$_{4}$ to be 1751 ppb, we solve for
$k_{M}$ as follows

$$k_{M} = \frac{-g\ln{(1-x)}}{0.554\cdot\nu\times10^{-9}P_{A}} = 1.4968$$

We then calculate $G_{M}$ as follows:

$$G_{M} = 0.554\times10^{-9}k_{M}\frac{P_{A}}{g} = 8.5410\times10^{-6}$$

Insolation
----------

From @mcgehee the latitude-dependence of insolation is given by the
relation,

$$insolation = Qs(y)$$

where $Q$ is the global annual average insolation and $y=sin(\varphi)$
with $\varphi$ representing latitude. $s(y)$ is a distribution of
relative insolation coefficients given by

$$s(y) = \frac{2}{\pi^{2}}\int_{0}^{2\pi}\sqrt{1-(
    \sqrt{1-y^{2}}\sin{\beta}\cos{\gamma}-y\cos{\beta})^{2}}d{\gamma}$$

$$\int_{0}^{1}s(y)dy = 1$$

where $\beta$ is the obliquity (here taken as $\beta$ = 23.5) and
$\gamma$ is the longitude. Taking the globally averaged solar radiation
as $Q$ = 340 W m$^{-2}$ and realizing that
$y = sin(\varphi) \implies  dy = cos(\varphi)d\varphi$, we compute the
average insolation over the Arctic region (A) defined as any region with
latitude $\geq$60as

$$\begin{aligned}
 Q_{A} & =  Q\cdot\frac{\int_{A}s(y)dy}{\int_{A}dy}\\
& =  Q\cdot \frac{\int_{\frac{\pi}{3}}^{\frac{\pi}{2}}\frac{2}{\pi^{2}}\int_{0}^{2\pi}\sqrt{1-(
    \sqrt{1-(\sin{\varphi})^{2}}\sin{\beta}\cos{\gamma}-\sin{\varphi}\cos{\beta})^{2}}d{\gamma}\cdot \cos(\varphi)d\varphi}{\int_{\frac{\pi}{3}}^{\frac{\pi}{2}}\cos(\varphi)d\varphi} \\
& \approx  201.73\:W\:m^{-2}\end{aligned}$$

