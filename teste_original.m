% for 20 time steps: 
t = 1:20; 
% starting locations of 30 points: 
x = 100*rand(30,1); 
y = 100*rand(30,1); 
% starting temperatures: 
T = 20+500*rand(30,1); 
% Make the first frame: 
figure
h = scatter(x+50*sind(t(1)/5),y+50*sind(t(1)/3),60,T*sind(t(1)/2),'filled');