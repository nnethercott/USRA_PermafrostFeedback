relinsol = @(x,y) 2.*sqrt(1-(sqrt(1-sin(y).^2).*sin(0.4101524).*cos(x)-sin(y).*cos(0.4101524)).^2)/pi.^2;
frac = @(y) integral(@(x) relinsol(x,y), 0, 2.*pi);

num = integral(@(y) frac(y)*cos(y), pi./3, pi./2, 'ArrayValued', true);
denom = integral(@(x) cos(x), pi./3, pi./2);

average_insol = 340*(num./denom)
%%
x=linspace(0,90,100);
for i=1:length(x)
    y(i) = frac(x(i).*pi./180);
end 

plot(x,y)
title('Insolation coefficient distribution')
xlabel('Latitude [deg N]')
ylabel('Coefficient')
ylim([0 1.4])
xlim([0 90])
yline(1, '--r')