% -------------------------------------------------------------------------
% Plot horizontal starting positions
% -------------------------------------------------------------------------

% Load starting points
[ start.lon start.lat start.p ] = textread('TEST','%f %f %f',-1); 

% Open a new figure and set the geographical projection and region
figure(1);
clf;
load coast
h=axesm('MapProjection','stereo','origin',[ 90 70 ]);
gridm;
h=plotm(lat,long,'Color','k','LineWidth',1.5)
%axis([-1.5 1.5 -2 0]);


% Plot starting points
for i=1:length(start.lon)
  linem(start.lat(i),start.lon(i),'marker','o', ...
        'markersize',4,'color','w','MarkerEdgeColor',[.4 .4 .4], ...
        'MarkerFaceColor','b');
end


% Save figure
figname = [ 'startf.eps' ];
set(gcf, 'PaperPosition', [2 1 15 10]);
print('-depsc2','-r0',figname);

