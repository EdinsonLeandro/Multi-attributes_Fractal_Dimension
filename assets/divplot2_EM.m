function [d,dn,m,R2] = divplot2_EM(X)
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

% delta=diff(X).^2 represent, for each row, the following equation:
% (a2-a1)^2 + (b2-b1)^2 + ... + (z2-z1)^2  (its depends on the number of
% coordinates in the data).
delta=diff(X).^2;

% It is "sum", because there is an "N" number of seismic attributes.

deltaTOTAL=sum(delta,2);
dmax=1/4*mean(deltaTOTAL);
dmin=sqrt(dmax);

fprintf('Computing -- please wait  ');

nn = 100;
i = 1;
d(i) = dmin;

% nn: number of steps for a given divisor. It changes from 2 (defined
% in by Middleton) to 3.
% In each cycle, divisor will increase to calculate the new steps number.

while nn > 2
    % Steps number to model the trace (depending on the divisor).
    nn = steps_EM(X,d(i));
    
    % Next divisor.
    dn(i) = nn*d(i);
    i = i + 1;
    d(i) = d(i-1)+dmin;
end;

% d = d(1:i-1,:);
d = d(1:i-1);
ld = log10(d);
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