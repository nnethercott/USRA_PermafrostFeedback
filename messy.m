g = ebm(9.75, 104, 201.73, 0.6);
g.mu = carbon(337);
g.nu = methane(337);

eqtemp(g,0.9)
g.tau_s
%%
fsolve(@(z) eqtempProfile(g,z)-1,0)