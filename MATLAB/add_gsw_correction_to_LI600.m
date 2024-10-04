function [data] = add_gsw_correction_to_LI600(filepath,stomatal_sidedness)
%ADD_GSW_CORRECTION_TO_LI600 Applies the Bailey & Rizzo (2024) correction
%of chamber air temperature and stomatal conductance to a csv file exported
%from an LI-600

% Input:
%  - filepath: Path to the CSV file exported from LI-600 (required).
%  - stomatal_sidedness: Correction factor for stomatal sidedness ...
%    1 if hypostomatous, 2 if amphistomatous, or anywhere in between (optional, default = 1).

% Output:
%  - new csv file with corrected gsw, T_chamber, W_chamber
if nargin < 2 || isempty(stomatal_sidedness)
    stomatal_sidedness = 1;  % Default
end

opts = detectImportOptions(filepath);
opts.VariableNamesLine = 2;
opts.DataLines = [4 inf];
data = readtable(filepath, opts);

sidedness = stomatal_sidedness*ones(size(data.gsw));
T_chambs = zeros(size(data.gsw));
T_outs = zeros(size(data.gsw));
W_chambs = zeros(size(data.gsw));
gsw_total = data.gsw.*sidedness;

for i=1:length(data.gsw)
    % --- inlet --- %
    T_in = data.Tref(i);                                % C
    RH_in = data.rh_r(i)/100;                           % Decimal 
    u_in = data.flow(i)/1000;                           % mmol/s
    P_atm = data.P_atm(i);                              % kPa
    
    % --- outlet --- %
    RH_out = data.rh_s(i)/100;                          % Decimal
    u_out = data.flow_s(i)/1000;                        % mmol/s, not used, deemed unreliable by LI-COR

    % --- chamber --- %
    T_leaf = data.Tleaf(i);                             % C
    RH_chamb = RH_out;                                  % Decimal

    % -- constants -- %
    s = 0.441786*0.01^2;                                % m^2
    gbw = 2.921;                                        % mol/m^2/s

    syms T_out;

    % -- calculate inlet values -- %
    W_in = calculateW(T_in,RH_in,P_atm);                % mol/mol
    h_in = calculateh(T_in,W_in);                       % J/mol
    Q_in = calculateQ(T_leaf,T_out,gbw);                % J/m^2/s
    
    % -- setup and solve implicit equation (4) from Bailey and Rizzo (2024) for T_out -- %
    W_out = calculateW(T_out,RH_out,P_atm);             % mol/mol
    h_out = calculateh(T_out,W_out);                    % J/mol 
    equation = (s.*Q_in == u_in*1000.*(h_out - h_in));
    T_chamb = vpasolve(equation, T_out);                % C
    T_out_value = T_chamb;

    % -- ASSUMPTION:The chamber air temperature is the average of the inlet and outlet air temperatures -- %
    T_chamb = 0.5*(T_in+T_chamb);                       % C
    
    W_chamb = calculateW(T_chamb,RH_chamb,P_atm);       % mol/mol
    W_leaf = calculateW(T_leaf,1,P_atm);                % mol/mol
    E = (u_in.*(W_chamb-W_in))./(s.*(1-W_chamb));       % mmol/m^2/s
    gtw = E./(W_leaf-W_chamb)/1000;                     % mol/m^2/s
    gsw_bottom(i) = 1./(1./gtw-1./gbw);                 % mol/m^2/s
    
    gsw_total(i) = gsw_bottom(i).*sidedness(i);         % mol/m^2/s
    T_chambs(i) = T_chamb;                              % C
    W_chambs(i) = W_chamb;                              % mol/mol
    T_outs(i) = T_out_value;                            % C
end

data.gsw_corrected = gsw_total;
data.Ta_chamb_corrected = T_chambs;
data.Wa_chamb_corrected = W_chambs;
data.T_out_corrected = T_outs;
data.stomatal_sidedness = sidedness;

[p,f,e]=fileparts(filepath);
filename=fullfile(p,f);
writetable(data,filename+"_corrected"+e);

end

function W = calculateW(T,RH,P)
a = 0.61365;
b = 17.502;
c = 240.97;
W = a*exp(b*T./(T+c)).*RH/P; % mol/mol
end

function h = calculateh(T,W)
cpa = 29.14; % J/mol/C
cpw = 33.5;  % J/mol/C
lambdaw = 45502; % J/mol
h = cpa*T+W.*(lambdaw + cpw*T);
end

function Q = calculateQ(T_in,T_out,gbw)
cpa = 29.14; % J/mol/C
heat_to_water = 1.08;
Q = cpa.*gbw./heat_to_water.*(T_in-T_out); % J/m^2/s
end