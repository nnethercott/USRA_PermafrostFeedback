r = rcpsim_8pool(10, 104, 201.73, 0.6);
r.ebm.Q_10 = 1.5;
%%
%interpolate to find max forcing 
maxvalFO = @(y) (y - 2.6).*(19.5 - 13)./(8.5 - 2.6) + 13;
maxvalFA = @(y) 129;
%maxvalFA = @(y) (y - 2.6).*(97.53 - 104.03)./(8.5 - 2.6) + 104.03;

%linear functions for each rcp 
linFO26 = @(x) 9.75 + (x-1900).*(maxvalFO(2.6)-9.75)./(2100-1900);
linFA26 = @(x) 107.28 + (x-1900).*(maxvalFA(2.6)-107.28)./(2100-1900);

linFO45 = @(x) 9.75 + (x-1900).*(maxvalFO(4.5)-9.75)./(2100-1900);
linFA45 = @(x) 107.28 + (x-1900).*(maxvalFA(4.5)-107.28)./(2100-1900);

linFO60 = @(x) 9.75 + (x-1900).*(maxvalFO(6.0)-9.75)./(2100-1900);
linFA60 = @(x) 107.28 + (x-1900).*(maxvalFA(6.0)-107.28)./(2100-1900);

linFO85 = @(x) 9.75 + (x-1900).*(maxvalFO(8.5)-9.75)./(2100-1900);
linFA85 = @(x) 107.28 + (x-1900).*(maxvalFA(8.5)-107.28)./(2100-1900);
times = linspace(1900, 2300, 401);

for i=1:(2100-1899)
    FO26(i) = linFO26(times(i));
    FA26(i) = linFA26(times(i));
    FO45(i) = linFO45(times(i));
    FA45(i) = linFA45(times(i));
    FO60(i) = linFO60(times(i));
    FA60(i) = linFA60(times(i));
    FO85(i) = linFO85(times(i));
    FA85(i) = linFA85(times(i));
end 

FO_arr26 = [FO26 maxvalFO(2.6).*ones([1 length(times)-length(FO26)])];
FA_arr26 = [FA26 maxvalFA(2.6).*ones([1 length(times)-length(FA26)])];
FO_arr45 = [FO45 maxvalFO(4.5).*ones([1 length(times)-length(FO26)])];
FA_arr45 = [FA45 maxvalFA(4.5).*ones([1 length(times)-length(FA26)])];
FO_arr60 = [FO60 maxvalFO(6.0).*ones([1 length(times)-length(FO26)])];
FA_arr60 = [FA60 maxvalFA(6.0).*ones([1 length(times)-length(FA26)])];
FO_arr85 = [FO85 maxvalFO(8.5).*ones([1 length(times)-length(FO26)])];
FA_arr85 = [FA85 maxvalFA(8.5).*ones([1 length(times)-length(FA26)])];
%% Constant approach 
FA_arr26 = 104.*ones([1 length(times)]);
FA_arr45 = 104.*ones([1 length(times)]);
FA_arr60 = 104.*ones([1 length(times)]);
FA_arr85 = 104.*ones([1 length(times)]);
FO_arr26 = 10.*ones([1 length(times)]);
FO_arr45 = 10.*ones([1 length(times)]);
FO_arr60 = 10.*ones([1 length(times)]);
FO_arr85 = 10.*ones([1 length(times)]);
%% GENERATE SURFACE TEMP CURVES 
%8.5
r.step = 0.25;
r.FO_arr = FO_arr85;
r.FA_arr = FA_arr85;
%t85e = eulerPermafrost(r, r.CO285data(1900-1764:2300-1764), r.CH485data(1900-1764:2300-1764));
info85p = rk3ab3(r, r.CO285data(1900-1764:2300-1764), r.CH485data(1900-1764:2300-1764));
t85 = no_permafrost(r, r.CO285data(1900-1764:2300-1764), r.CH485data(1900-1764:2300-1764));

%6.0
r.step = 0.5;
r.FO_arr = FO_arr60;
r.FA_arr = FA_arr60;
%t60e = eulerPermafrost(r, r.CO260data(1900-1764:2300-1764), r.CH460data(1900-1764:2300-1764));
info60p = rk3ab3(r, r.CO260data(1900-1764:2300-1764), r.CH460data(1900-1764:2300-1764));
t60 = no_permafrost(r, r.CO260data(1900-1764:2300-1764), r.CH460data(1900-1764:2300-1764));
%%
%4.5
r.FO_arr = FO_arr45;
r.FA_arr = FA_arr45;
%t45e = eulerPermafrost(r, r.CO245data(1900-1764:2300-1764), r.CH445data(1900-1764:2300-1764));
info45p = rk3ab3(r, r.CO245data(1900-1764:2300-1764), r.CH445data(1900-1764:2300-1764));
t45 = no_permafrost(r, r.CO245data(1900-1764:2300-1764), r.CH445data(1900-1764:2300-1764));

%2.6
r.FO_arr = FO_arr26;
r.FA_arr = FA_arr26;
%t26e = eulerPermafrost(r, r.CO226data(1900-1764:2300-1764), r.CH426data(1900-1764:2300-1764));
info26p = rk3ab3(r, r.CO226data(1900-1764:2300-1764), r.CH426data(1900-1764:2300-1764));
t26 = no_permafrost(r, r.CO226data(1900-1764:2300-1764), r.CH426data(1900-1764:2300-1764));
%%
t85p = downsample(info85p{19}, 4);
t60p = downsample(info60p{19}, 2);
t45p = downsample(info45p{19}, 2);
t26p = downsample(info26p{19}, 2);


t85 = t85(1:400);
t60 = t60(1:400);
t45 = t45(1:400);
t26 = t26(1:400);
%% PLOT 
times = linspace(1900, 2300, length(t85p));

%85
plot(times, t85p.*273.15 - 273.15, 'r--', 'LineWidth',1.2, 'DisplayName', 'RCP 8.5')
hold on
p1 = plot(times, t85.*273.15 - 273.15, 'b-', 'LineWidth',0.5, 'HandleVisibility', 'off');
p1.Color(4) = 0.6;

%60
plot(times, t60p.*273.15 - 273.15, 'g--', 'LineWidth',1.2, 'DisplayName', 'RCP 6.0')
p2 = plot(times, t60.*273.15 - 273.15, 'b-', 'LineWidth',0.5, 'HandleVisibility', 'off');
p2.Color(4) = 0.6;

%45
plot(times, t45p.*273.15 - 273.15, 'm--', 'LineWidth',1.2, 'DisplayName', 'RCP 4.5')
p3 = plot(times, t45.*273.15 - 273.15, 'b-', 'LineWidth',0.5, 'HandleVisibility', 'off');
p3.Color(4) = 0.6;

%26
plot(times, t26p.*273.15 - 273.15, 'k--', 'LineWidth',1.2, 'DisplayName', 'RCP 2.6')
p4 = plot(times, t26.*273.15 - 273.15, 'b-', 'LineWidth',0.5, 'DisplayName', 'No Permafrost');
p4.Color(4) = 0.6;

legend;
box off;
xlabel('Year')
ylabel('Arctic surface temperature [C]');


%% TESTING 

%8.5
r.FO_arr = FO_arr85;
r.FA_arr = FA_arr85;
info85 = rk3ab3INFO(r, r.CO285data(1900-1764:2300-1764), r.CH485data(1900-1764:2300-1764));

%6.0
r.FO_arr = FO_arr60;
r.FA_arr = FA_arr60;
info60 = rk3ab3INFO(r, r.CO260data(1900-1764:2300-1764), r.CH460data(1900-1764:2300-1764));

%4.5
r.FO_arr = FO_arr45;
r.FA_arr = FA_arr45;
info45 = rk3ab3INFO(r, r.CO245data(1900-1764:2300-1764), r.CH445data(1900-1764:2300-1764));

%2.6
r.FO_arr = FO_arr26;
r.FA_arr = FA_arr26;
info26 = rk3ab3INFO(r, r.CO226data(1900-1764:2300-1764), r.CH426data(1900-1764:2300-1764));
