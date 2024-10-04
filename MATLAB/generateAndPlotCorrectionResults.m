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