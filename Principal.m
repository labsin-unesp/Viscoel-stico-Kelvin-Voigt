%--------------------------------------------------------------------------   
%                 IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%               Universidade Estadual Paulista, FEIS/UNESP

%                  ESTUDO COMPORTAMENTAL DOS SOLOS
%        Autores: Proença, S. M.; Paschoalini, T. A; Campos, B. J.
%  Objetivo: Implementação do modelo Visco-Elástico Linear, com validação 
%      do domínio pela Função de Plastificação de Mohr-Coulomb
%                

%------------------ ANÁLISE TRANSIENTE (Excitação senoidal)----------------

np=3;                          % Numero de pontos de integracao (deve-se usar np>=3)
Nnos=(Ndivx+1)*(Ndivy+1);      % Calculando o número de nós
Nele=Ndivx*Ndivy;              % Calculando o número de elementos

% < GERAÇÃO AUTOMÁTICA DA MALHA >
disp([datestr(now,'<HH:MM:SS>') ' Iniciando geração da malha'])
[x,y,ice_N]=mesh(xv,yv,Ndivx,Ndivy,Nele);  

% Obs. Os nós são enumerados da esquerda para a direita de baixo para cima
% Ice_N = Matriz que possui em cada linha os nós que formam cada elemento (giro anti-horário)
% x e y são as coordenadas dos nós

figure(1)
for e=1:1:Nele
hold on
plot(x,y,'.b','MarkerSize',9);
xi=x(ice_N(e,1));
yi=y(ice_N(e,1));
xj=x(ice_N(e,2));
yj=y(ice_N(e,2));
xi2=x(ice_N(e,3));
yi2=y(ice_N(e,3));
xj2=x(ice_N(e,4));
yj2=y(ice_N(e,4));
k1=plot([xi xj],[yi yj],'k','linewidth',0.3);
k2=plot([xi2 xj2],[yi2 yj2],'k','linewidth',0.3);
k3=plot([xi xj2],[yj yj2],'k','linewidth',0.3);
k4=plot([xj xi2],[yj yi2],'k','linewidth',0.3);
end
title('Discretização da Geometria')
xlabel('X','FontSize',10)
ylabel('Y','FontSize',10)

[ice_G]=iceG(Nele,ice_N);      % Possui em cada linha os graus de liberdade de cada elemento (2gld por nó - "u" e "w")
[ide_L]=condcont(Ndivx,Ndivy); % Matriz formada pelos graus de liberdade livre

co1= 1000;                     % Coesão [Pa] = 1 kPa         
aatr1= 35;                     % Ângulo de atrito [°]
dila1= 9;                      % Ângulo de dilatância [°] (0.3*aatr)
%zeta=0.1;                      % Constante de amortecimento viscoso

E=E1*ones(1,Nele);             % Construindo o vetor de propriedades inicialmente com o valor da C1
v=v1*ones(1,Nele);      
co=co1*ones(1,Nele);    
aatr=aatr1*ones(1,Nele);    
dila=dila1*ones(1,Nele);     
p=p1*ones(1,Nele);     

h=20e-3*ones(1,Nele);          % Espessura para o EPT (obs: para o EPD assume valor unitário [m] - Direção mais longa)

%--------------------------------------------------------------------------
%               < INSERINDO O TERMO FONTE >


% ------------------ ENTRADA SENOIDAL --------------------
% Obs: Inicialmente, para simplificação, assumiremos uma excitação senoidal para o problema

rm=rms(input)                                 % Busca-se um rms próximo a 0.1
Finput=fft(input);
%input((floor(length(t)/3)):length(t)) = 0;   % Modificando o sinal de entrada (0<t<tf/2 - input=sin; tf/2<t<tf - 0)

% figure(2)
% subplot(2,1,1)
% plot(t,input)
% title('Excitação senoidal na base')
% xlabel('tempo de amostragem','FontSize',10)
% ylabel('Aceleração','FontSize',10)
% subplot(2,1,2)
% semilogy(linspace(0,Fs,length(abs(Finput*T))),abs(Finput*T),'r','linewidth',1)
% xlim([0 600])


%Neta=E1*neta;              % sidade dinâmica efetiva [Pa.s]   (7,6mPa.s < Neta < 500 kPa.s)  (usamos o E pq estamos "compac")

% Para o intervalo validado experimentalmente neta varia entre:
neta_min=7.6e-3/E1;
neta_max=500e3/E1;

% -------------------------------------------------------------------------
% < CRIANDO O VETOR DE FORÇA A PARTIR DA ACELERAÇÃO DE ENTRADA (TERMO FONTE) >
 
aci=zeros(2*Nnos,length(t));        % Cada coluna correponde a intensidade da aceleração num instante de tempo
%f=zeros(2*Nnos,length(t));         % Cada coluna correponde a intensidade da força num instante de tempo

% OBS: Se Ndivx for par o nó é o central, se for impar é o nó que antecede o centro
% OBS2!: Checar se realmente teremos no nó central esse rms desejado de aceleração (senão ajustar para)  

aci(local,:)=input(1,:);          % Nó central

PI=linspace(-1,1,2);              % Usado no segundo "for" do processamento (Def e Tensões)
%parte2=logical(0);                % Inclui ou ignora o calculo das deformações, tensões e ondas (true(1)/false(0))

%                            PROCESSAMENTO!
%**************************************************************************
disp([datestr(now,'<HH:MM:SS>') ' Iniciando o processamento'])
processamento                     % Rotina para o processamento dos dados  


% Obs: Para salvar as matrizes desejadas basta:
% save('matrizMb.txt','mg','-ascii') - salvar a matriz de massa
% save('matrizKb.txt','kg','-ascii') - salvar a matriz de rigidez
%deslo05=deslocamento(Gsup,:);
%save('neta05.txt','deslo05','-ascii')
