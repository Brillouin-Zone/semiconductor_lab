% Semiconductors
fname = 'C:\von_Server\ETH\BSc Physics\5\Praktikum 3\Semiconductor\plots';
%{
% 1. PART
%{
opts = spreadsheetImportOptions("NumVariables", 7);
opts.Sheet = "Tabelle1";
opts.DataRange = "A2:G23";
opts.VariableNames = ["TC", "ImuA", "V4mV", "err_V4mV", "V2mV", "err_V2mV", "err_I"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7], "EmptyFieldRule", "auto");
tbl = readtable("C:\von_Server\ETH\BSc Physics\5\Praktikum 3\Semiconductor\data\data.xlsx", opts, "UseExcel", false);
TC = tbl.TC;
ImuA = tbl.ImuA;
V4mV = tbl.V4mV;
err_V4mV = tbl.err_V4mV;
V2mV = tbl.V2mV;
err_V2mV = tbl.err_V2mV;
err_I = tbl.err_I;
clear opts tbl
    T = str2double(TC); % �C
    I = str2double(ImuA)*10^-6; % A
    V4 = str2double(V4mV)*10^-3; % V
    err_V4 = str2double(err_V4mV)*10^-3; % V
    V2 = str2double(V2mV)*10^-3; % V
    err_V2 = str2double(err_V2mV)*10^-3; % V
    err_T = 0.1*ones(22,1);% �C
    err_I = str2double(err_I)*10^-6; % A
    
    R_s = V4(:) ./ (I(:)); % range: (2:11)
    err_R_s = zeros(22,1);
    for i=1:22
        err_R_s(i) = sqrt((err_V4(i) / I(i))^2+(V4(i) * err_I(i) / I(i)^2)^2);
    end       
   
% CURRENT-LIMIT-PLOT    
%{
 figure
 yyaxis left
 errorbar(I(:), R_s(:), err_R_s(:)/2, err_R_s(:)/2, err_I(:)/2, err_I(:)/2, 'b*')
 ylabel('R_s [\Omega]');
 hold on
 h = vline([0.1*10^-3, 0.000555, 0.3*10^-3], {'k', 'k', 'k'}, {' ', ' ', ' '});
 hold on
 X = [0.0001, 0.000555, 0.0003];
 Y = [2815, 2815, 2815];
 str = {'\leftarrow I_{min}', '\leftarrow I_{max}', '\leftarrow I_{chosen}'};
 text(X, Y, str);
 hold on
 yyaxis right
 errorbar(I(:), T(:), err_T(:)/2, err_T(:)/2, err_I(:)/2, err_I(:)/2, 'r*')
 grid on
 title('current / temperature vs. internal resistance');
 xlabel('I [A]');
 ylabel('T [^\circ C]');
 saveas(gcf, fullfile(fname, 'Rs_vs_I.eps'), 'epsc');
 %}
 
%}
% 2. PART: HEAT UP

opts = spreadsheetImportOptions("NumVariables", 6);
opts.Sheet = "Tabelle1";
opts.DataRange = "A2:F117";
opts.VariableNames = ["T_C_up", "errT_up", "V4_up", "errV4_up", "I_mA_up", "errI_mA_up"];
opts.VariableTypes = ["double", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [2, 3, 4, 5, 6], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 3, 4, 5, 6], "EmptyFieldRule", "auto");
tbl = readtable("C:\von_Server\ETH\BSc Physics\5\Praktikum 3\Semiconductor\data\data0.xlsx", opts, "UseExcel", false);
    T_oven_up = tbl.T_C_up + 273.15; % K
    errToven_up = tbl.errT_up;
    V4oven_up = tbl.V4_up;
    errV4oven_up = tbl.errV4_up;
    I_mAoven_up = tbl.I_mA_up;
    errI_mAoven_up = tbl.errI_mA_up;
    clear opts tbl
    
    err_T_oven_up = str2double(errToven_up); % �C or K
    V4_oven_up = str2double(V4oven_up) * 10^-3; % V
    err_V4_oven_up = str2double(errV4oven_up)* 10^-3; % V
    I_oven_up = str2double(I_mAoven_up)* 10^-3; % A
    err_I_oven_up = str2double(errI_mAoven_up)* 10^-3; % A
    
    R_s_oven_up = V4_oven_up(:) ./ I_oven_up(:);
    err_R_s_oven_up = zeros(99,1);
    for i= 1:99
        err_R_s_oven_up(i) = sqrt((err_V4_oven_up(i) / I_oven_up(i))^2+(V4_oven_up(i) * err_I_oven_up(i) / I_oven_up(i)^2)^2);
    end
         
    figure
    %errorbar(T_oven_up(:), R_s_oven_up(:), err_R_s_oven_up(:)/2, err_R_s_oven_up(:)/2, err_T_oven_up(:)/2, err_T_oven_up(:)/2, 'b*')
    plot(T_oven_up(:), R_s_oven_up(:), 'b*')
    X1 = [280, 430, 650];
    Y1 = [3500, 2500, 300];
    str1 = {'1', '2', '3'};
    text(X1, Y1, str1)
    xlabel('temperature T [K]');
    ylabel('resistance R_s [\Omega]');
    title('temperature-dependent resistance');
    saveas(gcf, fullfile(fname, 'oven_up.eps'), 'epsc');
   
    lnR_up = log(R_s_oven_up(:));
    kB = 8.617330350*10^-5; %eV/K
    X_up = (2*kB * T_oven_up(:)).^-1;
        % from curvefit.m
        %[xXData, yYData] = prepareCurveData( X_up, lnR_up ); 
        [xXData, yYData] = prepareCurveData( X_up(16:116), lnR_up(16:116) ); % ohne rechteste 15 data points. wegen 1/T sind dies die ersten 15
        fFt = fittype( 'poly1' );
        [fFitresult, gGof] = fit( xXData, yYData, fFt );
            % fFitresult: p1 =       1.036  (1.031, 1.04)
            %               -> err_E_gap_heating := 0.005
            %             p2 =       -8.04  (-8.085, -7.995)
        cCoeff = coeffvalues(fFitresult);
        fFitt1 = fit(xXData, yYData, fFt);
        
                
    figure
    %plot( fFitresult, xXData, yYData);
    %hold on
    %errorbar(X(:), lnR(:), err_R_s_oven(:)/2, err_R_s_oven(:)/2, err_T_oven(:)/2, err_T_oven(:)/2, 'b*')
    plot(X_up(:), lnR_up(:), 'b.')
    hold on
    plot(X_up(16:116), cCoeff(1)*X_up(16:116) + cCoeff(2), 'r-')
    hold on
    xX = 8.566;
    yY = 0.8866;% (x,y) beschreibt den punkt, wo die Fitline angeschrieben werden soll
    sStr = {'\leftarrow y = 1.0358x -8.0400'};
    text(xX, yY, sStr, 'Color', 'red');
    X3 = [11, 18];
    Y3 = [4, 7.5];
    str3 = {'1', '2'};
    text(X3, Y3, str3)
    xlabel('(2k_B T)^{-1} [a.u.]');
    ylabel('ln(R_s) [a.u.]');
    legend('data', 'linear fit', 'Location', 'Best');
    title('T-dependent resistance while heating up');
    saveas(gcf, fullfile(fname, 'oven_lnR_up.eps'), 'epsc');
    
  %{
            % Figure to fit linearly for residuals
            figure
            plot(X_up(16:116), lnR_up(16:116), 'bx')
            hold on
            plot(X_up(16:116), cCoeff(1)*X_up(16:116) + cCoeff(2), 'r-')
            hold on
            xX = 8.566;
            yY = 0.8866;% (x,y) beschreibt den punkt, wo die Fitline angeschrieben werden soll
            sStr = {'\leftarrow y = 1.0358x -8.0400'};
            text(xX, yY, sStr, 'Color', 'red');
            X3 = [11, 18];
            Y3 = [4, 7.5];
            str3 = {'1', '2'};
            text(X3, Y3, str3)
            xlabel('(2k_B T)^{-1} [a.u.]');
            ylabel('ln(R_s) [a.u.]');
            legend('data', 'linear fit', 'Location', 'Best');
            title('FOR RESIDUALS heating');
            saveas(gcf, fullfile(fname, 'oven_lnR_up_FOR_RESIDUALS.eps'), 'epsc');
    
    %}
E_gap_up = 2*kB * T_oven_up(:) .* lnR_up(:);
    
    
% 2. PART : COOL DOWN

opts = spreadsheetImportOptions("NumVariables", 6);
opts.Sheet = "Tabelle1";
opts.DataRange = "A2:F100";
opts.VariableNames = ["T_C", "errT", "V4", "errV4", "I_mA", "errI_mA"];
opts.VariableTypes = ["double", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [2, 3, 4, 5, 6], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 3, 4, 5, 6], "EmptyFieldRule", "auto");
tbl = readtable("C:\von_Server\ETH\BSc Physics\5\Praktikum 3\Semiconductor\data\data2.xlsx", opts, "UseExcel", false);
    T_oven = tbl.T_C + 273.15; % K
    errToven = tbl.errT;
    V4oven = tbl.V4;
    errV4oven = tbl.errV4;
    I_mAoven = tbl.I_mA;
    errI_mAoven = tbl.errI_mA;
    clear opts tbl
    
    err_T_oven = str2double(errToven); % �C or K
    V4_oven = str2double(V4oven) * 10^-3; % V
    err_V4_oven = str2double(errV4oven)* 10^-3; % V
    I_oven = str2double(I_mAoven)* 10^-3; % A
    err_I_oven = str2double(errI_mAoven)* 10^-3; % A
    
    R_s_oven = V4_oven(:) ./ I_oven(:);
    err_R_s_oven = zeros(99,1);
    for i= 1:99
        err_R_s_oven(i) = sqrt((err_V4_oven(i) / I_oven(i))^2+(V4_oven(i) * err_I_oven(i) / I_oven(i)^2)^2);
    end
        
    figure
    %errorbar(T_oven(:), R_s_oven(:), err_R_s_oven(:)/2, err_R_s_oven(:)/2, err_T_oven(:)/2, err_T_oven(:)/2, 'b*')
    plot(T_oven(:), R_s_oven(:), 'b*')
    X2 = [430, 650];
    Y2 = [2500, 300];
    str2 = {'2', '3'};
    text(X2, Y2, str2)
    xlabel('temperature T [K]');
    ylabel('resistance R_s [\Omega]');
    title('temperature-dependent resistance');
    saveas(gcf, fullfile(fname, 'oven.eps'), 'epsc');
    
    lnR = log(R_s_oven(:));
    kB = 8.617330350*10^-5; %eV/K
    X = (2*kB * T_oven(:)).^-1;
        % from curvefit.m
        [xData, yData] = prepareCurveData( X, lnR );
        ft = fittype( 'poly1' );
        [fitresult, gof] = fit( xData, yData, ft );
            % fitresult: p1 = 1.113  (1.107, 1.118) 
            %            -> err_E_gap_cooling := 0.005
            %            p2 = -8.647  (-8.704, -8.589)
            
        coeff = coeffvalues(fitresult);
        fitt1 = fit(xData, yData, ft);

    figure
    plot( fitresult, xData, yData);
    hold on
    %plot(X(:), lnR(:), 'bo')
    %errorbar(X(:), lnR(:), err_R_s_oven(:)/2, err_R_s_oven(:)/2, err_T_oven(:)/2, err_T_oven(:)/2, 'b*')
    x = 8.566; %
    y = 0.8866; % (x,y) beschreibt den punkt, wo die Fitline angeschrieben werden soll
    str = {'\leftarrow y = 1.113x -8.646'};
    text(x, y, str, 'Color', 'red');
    X4 = [11, 15.5];
    Y4 = [4, 7.8];
    str4 = {'1', '2'};
    text(X4, Y4, str4)
    xlabel('(2k_B T)^{-1} [a.u.]');
    ylabel('ln(R_s) [a.u.]');
    legend('data', 'linear fit', 'Location', 'Best');
    title('T-dependent resistance while cooling down');
    saveas(gcf, fullfile(fname, 'oven_lnR.eps'), 'epsc');
%}
%{        
    % for residuals
         figure
         plot(X(:), lnR(:), 'bo')
         xlabel('(2k_B T)^{-1} [a.u.]');
         ylabel('ln(R_s) [a.u.]');
         legend('data', 'linear fit', 'Location', 'Best');
         title('FOR RESIDUALS COOLING');
         saveas(gcf, fullfile(fname, 'oven_lnR_RESIDUALS_COOLING.eps'), 'epsc');
    
  %}  
  %}  
  
  
% CARRIER DENSITY AND MOBILITY

% Mobility: mu propto T^{-3/2}
mu_up = (T_oven_up(:)).^(-3/2);
err_mu_up = (3/2) *(T_oven_up(:)).^(-5/2) .* err_T_oven_up(:);
    figure
    errorbar(T_oven_up(:), mu_up(:), err_mu_up(:)/2, err_mu_up(:)/2, err_T_oven_up(:)/2, err_T_oven_up(:)/2, 'b-');
    %hold on
    %plot(T_oven_up(:), (T_oven_up(:)).^(-3/2), 'r-')
    xlabel('T [K]');
    ylabel('\mu [m^2 / Vs]');
    %legend('data', 'T^{-3/2}', 'Location', 'Best');
    title('Mobility \mu as a function of T');
    saveas(gcf, fullfile(fname, 'mobility.eps'), 'epsc');
    
 % carrier concentration n
 E_gap_heat = 1.04 * 1.602*10^-19; %V
 err_E_gap_heat = 0.005* 1.602*10^-19; %V
 E_gap_cool = 1.11 * 1.602*10^-19; %V
 err_E_gap_cool = 0.005* 1.602*10^-19; %V
 n_heat = ((T_oven_up(:)).^(3/2)).' .* exp(-E_gap_heat / (2*kB .* T_oven_up(:)));
 err_n_heat = zeros(1,116);
 n_cool = ((T_oven(:)).^(3/2)).' .* exp(-E_gap_cool / (2*kB .* T_oven(:)));
 err_n_cool = zeros(1,99);
 for i=1:116
     err_n_heat(i) = ((err_T_oven_up(i) * (3/2 * (T_oven_up(i))^(1/2)*exp(-E_gap_heat / (2*kB * T_oven_up(i)))+ (T_oven_up(i))^(3/2) * (E_gap_heat / 2*kB *(T_oven_up(i))^2)*exp(-E_gap_heat / (2*kB * T_oven_up(i)))))^2 + (err_E_gap_heat * (T_oven_up(i))^(3/2) * exp(-E_gap_heat / (2*kB * T_oven_up(i))) * (1/(2*kB*T_oven_up(i))))^2)^(1/2);
 end
 for j=1:99
    err_n_cool(j) = ((err_T_oven(j) * (3/2 * (T_oven(j))^(1/2)*exp(-E_gap_cool / (2*kB * T_oven(j)))+ (T_oven(j))^(3/2) * (E_gap_cool / 2*kB *(T_oven(j))^2)*exp(-E_gap_cool / (2*kB * T_oven(j)))))^2 + (err_E_gap_cool * (T_oven(j))^(3/2) * exp(-E_gap_cool / (2*kB * T_oven(j))) * (1/(2*kB*T_oven(j))))^2)^(1/2);
 end
 
 
    figure
    %plot(T_oven_up(:), n_heat(:), 'b-', T_oven(:), n_cool(:), 'c-');
    %plot(T_oven_up(:), n_heat(:), 'b-');
    errorbar(T_oven_up(:), n_heat(:), err_n_heat(:)/2, err_n_heat(:)/2, err_T_oven_up(:)/2, err_T_oven_up(:)/2, 'b-');
    %legend('heating process', 'cooling process');
    title('Charge carrier density n as a function of T');
    xlabel('T [K]');
    ylabel('n [1 / m^3]');
    saveas(gcf, fullfile(fname, 'carrier_density.eps'), 'epsc');
    
   n_Si_300K = 300^(3/2) * exp(-E_gap_heat / (2*kB * 300))
