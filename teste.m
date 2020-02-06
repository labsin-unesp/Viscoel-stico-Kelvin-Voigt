% for 20 time steps: 
t = 1:20; 
% starting locations of 30 points: 

%x = 100*rand(30,1); 
%y = 100*rand(30,1);
[x,y,ice_N]=mesh(xv,yv,Ndivx,Ndivy,Nele);



% starting temperatures: Primeira coluna da matriz
T = campo_t(Ndivx,Ndivy,aceleracao,1);

% Make the first frame: 
figure
h = scatter(x+50*sind(t(1)/5),y+50*sind(t(1)/3),60,T*sind(t(1)/2),'filled'); 
% set x/y axis limits: 
axis([0 100 0 100]) 
% set color axis limits: 
caxis([20 25]) 
cb = colorbar; 
ylabel(cb,'temperature')

% write the first frame: 
% gif('temperaturedata.gif') <-uncomment to make a gif
% Loop through each subsequent time step: 
for k = 2:length(t) 
     set(h,'xdata',x+50*sind(t(k)/5),'ydata',y+50*sind(t(k)/3),...
        'cdata',T*sind(t(k)/2))
     pause(0.1) % <-not necessary for making a gif
     %gif%  <-uncomment to make a gif
 end