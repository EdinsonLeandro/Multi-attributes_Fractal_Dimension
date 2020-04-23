function graph_QC(Traces,nPD_R2,rL,COORD,TIME,Title,i)
% Make QC plots of Fractal Dimension calculations.
% Traces = Seismic Attributes Traces selected from the time window of
%          interest according to the surface.
% nPD_R2 = Matrix: Selected trace number / Slope of the line
%          between log10(r) and log10(L) / Fractal Dimension / R2.
%     rL = Hypermatrix. Row 1: Dividers. Row 2: Total Length. For Fractal
%          Dimension calculations. The 3rd Dimension is each trace.
%  COORD = Matrix with the coordinates of each trace selected.
%   TIME = Matrix with time windows for each D calculation.
%  Title = Cell array with the X-Y labels of seismic attributes 3D graph.
%      i = Input trace index for QC.

figure(i)

% 3D graph of two traces. If two or more seismic attributes were used in 
% the Fractal Dimension calculation, only the first two will be used here.

trace_number = num2str(nPD_R2(i,1));
x = num2str(COORD(i,1));
y = num2str(COORD(i,2));
if size(Traces,2)==1
    xytitle = ['Traza de Atributo Sísmico #' trace_number ...
        '.  Coordenadas =  E: ' x '   N: ' y];
    subplot(1,2,1);
    plot(Traces(:,1,i), TIME(:,i),'b','LineWidth',2);
    title(xytitle);
    grid
    set(gca,'YDir','reverse', 'LineWidth',2,'Color','w')
    axis('auto')
    xlabel(Title{1})
    ylabel('Tiempo')
else
    xytitle = ['Vista 3D Atributos Sísmicos. Traza #' trace_number ...
        '.  Coordenadas =  E: ' x '   N: ' y];
    subplot(1,2,1);
    plot3(Traces(:,1,i), Traces(:,2,i), TIME(:,i),'b','LineWidth',2);
    title(xytitle);
    grid
    set(gca,'ZDir','reverse', 'LineWidth',2,'Color','w')
    axis('vis3d')
    xlabel(Title{1})
    ylabel(Title{2})
    zlabel('Tiempo')
end

% Plot: Fractal Dimension.

slope = num2str(nPD_R2(i,2));
ss = ['Traza #' trace_number '     Pendiente = ' slope];

% Logarithm with base 10, for divisor and total length. As a previous step,
% it must eliminate the possible content of zeros within the matrix
% (which are generated due to variability of matrix dimensions).

r=rL(1,:,i);
L=rL(2,:,i);
empty_rL = find(r==0);
r(empty_rL) = [];
L(empty_rL) = [];
ld = log10(r);
ldn = log10(L);

% Linear regression
p = polyfit(ld,ldn,1);
f = polyval(p,ld);

hold on
subplot(1,2,2);
plot(ld,ldn,'ks',ld,f,'g','LineWidth',2,'MarkerFaceColor','b');
axis('equal');
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
text(mean(xlim), 0.7*ylim(2),['R^2 = ', num2str(nPD_R2(i,4))],...
    'FontWeight','bold','HorizontalAlignment','center');
ylabel('Log10 (Longitud Total)'), xlabel('Log10 (Divisor)');
title(ss);


