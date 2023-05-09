function [dividers,dn,m,R2] = divplot2_EM(X)
% The original function: [d,dn] = divplot2(X)
%
% plots a divider plot, log(dn) vs log(d) for divider lengths, d,
% that are powers of 2 of dmin. Where:
% * dmin is twice the maximum distance between succesive points.
% * dn is the total measured distance for each divider length d.
% 
% Calls steps3.m 
% This version divplot2 uses a random row in the matrix X as starting point, 
% and computes dn both above and below this row: the sum of the two is the
% value returned. This is repeated for 10 different starting points.
% 
% Written by Gerry Middleton, McMaster University, July 1996
% EMail: middleto@mcmaster.ca

% This is a modified version:
% Distance calculation between each points.
% Attention! Two definitions of the divisors are setting:
% Middelton: dmin = 2*sqrt(dmax).
% Nellyana:  dmin=(1/8)*sqrt(dmean)

% dmin=minDivider.  dmax=maxDivider.  d=dividers.

% diff(X), for a vector X, is [X(2)-X(1), X(3)-X(2), ..., X(n)-X(n-1)]. 
% diff(X), for a matrix X, is the same calculation along each column.

% Therefore:
% diff(X).^2 = [(a(2)-a(1)).^2, (a(3)-a(2)).^2, ..., (a(n)-a(n-1)).^2;
%               (b(2)-b(1)).^2, (b(3)-b(2)).^2, ..., (b(n)-b(n-1)).^2; ...
%               (z(2)-z(1)).^2, (z(3)-z(2)).^2, ..., (z(n)-z(n-1)).^2]
% Where "a, b, ..., z" are the coordinates of the data.
delta = diff(X).^2;

% "sumDelta" represent, for each row, the following equation:
% (a2-a1)^2 + (b2-b1)^2 + ... + (z2-z1)^2
sumDelta = sum(delta,2);

% Maximun and minimun divider.
maxDivider = 1/4*mean(sumDelta);
minDivider = sqrt(maxDivider);

i = 1;

% Vector with divider used in each loop.

dividers(i) = minDivider;
nSteps = Inf;

% nSteps: Number of steps to approximate data depending on the given divisor.
% It changes from 2 (defined in by Middleton) to 3.


% In each cycle, divisor will increase to calculate the new steps number.

while nSteps > 2
    % Number of steps.
    nSteps = steps_EM(X, dividers(i));
    
    % Next divisor.
    dn(i) = nSteps*dividers(i);
    i = i + 1;
    dividers(i) = dividers(i-1) + minDivider;
end

% dividers = dividers(1:i-1,:);
dividers = dividers(1:i-1);
ld = log10(dividers);
ldn = log10(dn);

p = polyfit(ld,ldn,1); % linear regression
f = polyval(p,ld);
m = p(1);   % The slope is associated with Fractal Dimension(D) through
            % the equation: D=1+abs(m).

% R^2 calculation.
R2 = rsquare(ld,ldn,p);
%            
% % Plot.
% s = num2str(p(1));
% ss = ['Slope = ' s];
% plot(ld,ldn,'o',ld,f), axis('equal');
% 
% % Write R2 at the top of the graph.
% xlim = get(gca,'XLim');
% ylim = get(gca,'YLim');
% text('Position',[mean(xlim) 0.85*ylim(2)],'FontWeight','bold',...
%     'String',['R2 = ', num2str(R2)], 'HorizontalAlignment','center')
% 
% ylabel('log2(total length)'), xlabel('log2(interval)')
% title(ss);