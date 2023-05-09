function [R2] = rsquare(x, y, p)
% Coefficient of Determination in a linear regression model.
% Only Ordinary (unadjusted) R-squared.
% x = Independent Variable.
% y = Dependent Variable.
% p = Polynomial curve fitting coefficients, following the equation:
% p(x)= p(1)*x^n + p(2)*x^(n-1) + ... + p(n)*x + p(n+1).
%
% n=1 in this case. Therefore: p = [p1 , p2]
%
% Read more:
% www.mathworks.com/help/stats/coefficient-of-determination-r-squared.html

% The Matlab function "fitlm", creates a linear regression model more
% complete, with the calculation of all variables, including R^2.
% However, for efficiency reasons, it will not be used.

yObserved = y;
yTheoretical = (p(1).*x) + p(2);
meanY = mean(yObserved);

% SSE is the sum of squared error.
% SST is the sum of squared total.
SSE = sum((yObserved - yTheoretical).^2);
SST = sum((yObserved-meanY).^2);

R2 = 1 - (SSE/SST);