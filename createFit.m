function [fitresult, gof] = createFit(X, lnR)
%CREATEFIT(X,LNR)
%  Create a fit.
%
%  Data for 'lnR_fit' fit:
%      X Input : X
%      Y Output: lnR
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 20-Nov-2018 20:06:48


%% Fit: 'lnR_fit'.
[xData, yData] = prepareCurveData( X, lnR );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'lnR_fit' );
h = plot( fitresult, xData, yData );
legend( h, 'lnR vs. X', 'lnR_fit', 'Location', 'NorthEast' );
% Label axes
xlabel X
ylabel lnR
grid on

