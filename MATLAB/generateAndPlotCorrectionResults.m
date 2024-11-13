%% Generate Correction
clear;
InputFile = "walnut.csv";
data = add_gsw_correction_to_LI600(InputFile);
%% Plot
figure();
x = data.gsw;
y = data.gsw_corrected;
scatter(x,x,"k"); hold on;
scatter(x,y,"k","filled");
[res,gof] = fit(x,y,"a*x+c");
plot(res);
legend("Original g$_{sw}$","Corrected g$_{sw}$","y = "+round(res.a,2)+"x+"+round(res.c,4),"location","southeast","Interpreter","latex");
xlabel("Original g$_{sw}$ (mol m$^{-2}$ s$^{-1}$)","Interpreter","latex");
ylabel("g$_{sw}$ (mol m$^{-2}$ s$^{-1}$)","Interpreter","latex");
title("LI-600 Stomatal Correction","Interpreter","latex");
set(gca,"LineWidth",2);
set(gca,"TickLabelInterpreter","latex");
set(gca,"Color","white");
set(gca,"Fontsize",16)

%%
figure();
a = 0.61365;                                % unitless (empirical magnitude of es vs T)
b = 17.502;                                 % unitless (empirical slope of es vs T)
c = 240.97;                                 % C (empirical offset of es vs T)
es = @(T) a*exp(b*T./(T+c));                % kPa (saturation vapor pressure vs T function)
W = @(T,RH) es(T).*RH./(data.P_atm);        % mol/mol (water vapor mole fraction)
x = W(data.Tref,data.rh_r/100);
y = data.W_chamb_corrected;
scatter(x,x,"k"); hold on;
scatter(x,y,"k","filled");
[res,gof] = fit(x,y,"a*x+c");
plot(res);
legend("Original W$_{chamb}$","Corrected W$_{chamb}$","y = "+round(res.a,2)+"x+"+round(res.c,4),"location","southeast","Interpreter","latex");
xlabel("Original W$_{chamb}$ (mol m$^{-2}$ s$^{-1}$)","Interpreter","latex");
ylabel("W$_{chamb}$ (mol mol$^{-}$)","Interpreter","latex");
title("LI-600 Stomatal Correction","Interpreter","latex");
set(gca,"LineWidth",2);
set(gca,"TickLabelInterpreter","latex");
set(gca,"Color","white");
set(gca,"Fontsize",16)