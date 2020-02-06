%--------------------------------------------------------------------------   
%                 IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%              Universidade Estadual Paulista, FEIS/UNESP

%   Neste sub-algoritmo é realizado o processamento dos dados e a construção
% das matrizes da equação dinâmica: [M].{a} + [C].{v} + [K].{x} = [M].[I].{aex}

%---------------------------- PROCESSAMENTO -------------------------------

[PG,W]=quadgauss(np);        % Quadratura de gauss - Método numérico para resolução das integrais
u=zeros(2*Nnos,length(t));   % u = campo de deslocamentos
ve=zeros(2*Nnos,length(t));  % ve = campo de velocidades
ac=zeros(2*Nnos,length(t));  % ac = campo de acelerações
kg_xx=zeros(2*Nnos);         % kg_xx = matriz de rigidez global referente a xx
kg_yy=zeros(2*Nnos);         % kg_yy = matriz de rigidez global referente a yy
kg_xy=zeros(2*Nnos);         % kg_xy = matriz de rigidez global referente a xy
mg=zeros(2*Nnos);            % mg = matriz de massa global
fvg=zeros(Nnos,1);           % fvg = matriz global referente a ação gravitacional
gamag=zeros(Nnos);         % Parcela referente aos elementos em y=0 ou y=H

kg_yyOK=zeros(Nnos);kg_xxOK=zeros(Nnos); kg_xyOK=zeros(Nnos); mgOK=zeros(Nnos);

% CONSTRUINDO AS MATRIZES DO PROBLEMA
disp([datestr(now,'<HH:MM:SS>') ' Construindo as matrizes do problema'])
for e=1:1:Nele
     x11=x(ice_N(e,1));      % Coordenada x do nó 1
     x21=x(ice_N(e,2));      % Coordenada x do nó 2
     x31=x(ice_N(e,3));      % Coordenada x do nó 3
     x41=x(ice_N(e,4));      % Coordenada x do nó 4
     x12=y(ice_N(e,1));      % Coordenada y do nó 1
     x22=y(ice_N(e,2));      % Coordenada y do nó 2
     x32=y(ice_N(e,3));      % Coordenada y do nó 3
     x42=y(ice_N(e,4));      % Coordenada y do nó 4
     [J11i,J12i,J21i,J22i,detJ]=jacobiano(x11,x21,x31,x41,x12,x22,x32,x42,PG);
     [B11,B12,B13,B14,B15,B16,B17,B18,B21,B22,B23,B24,B25,B26,B27,B28,B31,B32,B33,B34,B35,B36,B37,B38]=matrizB(J11i,J12i,J21i,J22i,PG);
     k_xx=matrizK_xx(h(e),detJ,B11,B12,B13,B14,B15,B16,B17,B18,np);
     %kg_xx(ice_G(e,:),ice_G(e,:))=kg_xx(ice_G(e,:),ice_G(e,:))+k_xx;
     k_yy=matrizK_yy(h(e),detJ,B21,B22,B23,B24,B25,B26,B27,B28,np);
     %kg_yy(ice_G(e,:),ice_G(e,:))=kg_yy(ice_G(e,:),ice_G(e,:))+k_yy;
     k_xy=matrizK_xy(h(e),detJ,B11,B12,B13,B14,B15,B16,B17,B18,B21,B22,B23,B24,B25,B26,B27,B28,np);
     %kg_xy(ice_G(e,:),ice_G(e,:))=kg_xy(ice_G(e,:),ice_G(e,:))+k_xy;
     [N11,N12,N13,N14,N15,N16,N17,N18,N21,N22,N23,N24,N25,N26,N27,N28]=matrizN(np);
     m=massa(N11,N13,N15,N17,N22,N24,N26,N28,p(e),h(e),detJ,np);
     mg(ice_G(e,:),ice_G(e,:))=mg(ice_G(e,:),ice_G(e,:))+m;
     gama=contornos(h(e),detJ,N11,N12,N13,N14,N15,N16,N17,N18,B21,B22,B23,B24,B25,B26,B27,B28,np);
     gamag(ice_N(e,:),ice_N(e,:))=gamag(ice_N(e,:),ice_N(e,:))+gama((mod(ice_G(e,:),2)==1),(mod(ice_G(e,:),2)==0));
     
     kg_xxOK(ice_N(e,:),ice_N(e,:))=kg_xxOK(ice_N(e,:),ice_N(e,:))+k_xx((mod(ice_G(e,:),2)==1),(mod(ice_G(e,:),2)==1));
     kg_yyOK(ice_N(e,:),ice_N(e,:))=kg_yyOK(ice_N(e,:),ice_N(e,:))+k_yy((mod(ice_G(e,:),2)==0),(mod(ice_G(e,:),2)==0));
     kg_xyOK(ice_N(e,:),ice_N(e,:))=kg_xyOK(ice_N(e,:),ice_N(e,:))+k_xy((mod(ice_G(e,:),2)==0),(mod(ice_G(e,:),2)==1));
     mgOK(ice_N(e,:),ice_N(e,:))=mgOK(ice_N(e,:),ice_N(e,:))+m((mod(ice_G(e,:),2)==1),(mod(ice_G(e,:),2)==1));
     fv=gravidade(N11,N13,N15,N17,p,h,detJ,np);
     fvg(ice_N(e,:),1)=fvg(ice_N(e,:),1)+fv;
    
     % spy(kg_xxOK) % mostra a disposição dos elementos não nulos
end

% Após obter todas as matrizes globais do problema é possível solucionar a equação dinâmica:

% [coef11]{v} + [coef12]{w} = f1
% [coef21]{v} + [coef22]{w} = f2

elem_cont=[(1:Ndivx),((Nele-Ndivx+1):Nele)];   % Elementos que se encontram em y=0 ou y=H (fronteiras)
ide_Lno=ide_L(mod(ide_L,2)==0)./2;             % Nós que estão livres 
delt=dt;                                       % Delta t

% As matrizes multiplitativas (coef) são constantes no tempo

const1=((E1/(1-v1^2)).*kg_xxOK+G.*kg_yyOK);
const2=(((E1*v1/(1-v1^2))+G).*kg_xyOK);
const3=(G.*kg_xxOK+(E1/(1-v1^2)).*kg_yyOK);

coef11=((4/delt^2).*mgOK+((2*neta/delt)+1).*const1);
coef12=(((2*neta/delt)+1).*const2);
coef21=coef12;
coef22=((4/delt^2).*mgOK+((2*neta/delt)+1).*const3);

% -------------------------------------------------------------------------
% FORMATANDO AS MATRIZES:

% Construindo a Matriz A que será utilizada no método numérico de solução:
% [A]*{deslocamento}={F}

% onde: [A] = [[coef11] [coef12]] ; resul={v}  ; F={f1}
%              [coef21] [coef22]          {w}      {f2}

for linha=1:Nnos
  for coluna=1:Nnos
     sub_m=[coef11(linha,coluna) coef12(linha,coluna); coef21(linha,coluna) coef22(linha,coluna)];  % sub-matriz (2x2)
     mcell_A{linha,coluna}=sub_m;  % controi-se uma matriz celular com as submatrizes
  end
end    
matriz_A=cell2mat(mcell_A);

%--------------------------------------------------------------------------
% Precondicionador: (Acelera a convergência do método iterativo!)
[Ls,Us] = ilu(sparse(matriz_A(ide_L,ide_L)),struct('type','ilutp','droptol',1e-6));     % Fatoração Incompleta LU (ILU)
L=full(Ls);
U=full(Us);

% -------------------------------------------------------------------------
% MÉTODO DE SOLUÇÃO NUMÉRICO  - GMRES
disp([datestr(now,'<HH:MM:SS>') ' Iniciando método de solução'])

contt=1;                                       % contagem dos steps 

% Iniciando os vetores de deslocamento, velocidade e aceleração para o tempo anterior (n)
de_vo=zeros(Nnos,1); de_wo=zeros(Nnos,1);      % Vetores dos deslocamentos (tempo anterior)
ve_vo=zeros(Nnos,1); ve_wo=zeros(Nnos,1);      % Vetores das velocidades (tempo anterior)
ac_vo=zeros(Nnos,1); ac_wo=zeros(Nnos,1);      % Vetor das acelerações (tempo anterior)

ve_x=zeros(Nnos,length(t)); ve_y=zeros(Nnos,length(t)); % "Vetores" onde se armazenará as velocidades calculadas para cada instante
ac_x=zeros(Nnos,length(t)); ac_y=zeros(Nnos,length(t)); % "Vetores" onde se armazenará as acelerações calculadas para cada instante

% Vetor final com as duas váriaveis imbutidas e alternadas (v1; w1; v2; w2; ...)
deslocamento=zeros(2*Nnos,length(t));

converg=zeros(length(1:length(t)),1);

Gsup=2*(((Ndivx+1)*Ndivy))+(local); % Grau de liberdade do nó superficial acima do vazamento (em y)
ini = 1;
% Laço para solução do problema em cada instante de tempo
for inst=1:length(t)
   
   % f1, f2, e, consequentemente, F devem ser atualizados a cada passo de tempo!
   f1=(((4/delt^2).*mgOK+(2*neta/delt).*(const1))*de_vo)+((4/delt).*mgOK+neta.*const1)*ve_vo+mgOK*ac_vo+((2*neta/delt).*const2)*de_wo+neta.*const2*ve_wo;  
   f2=(((4/delt^2).*mgOK+(2*neta/delt).*(const3))*de_wo)+((4/delt).*mgOK+neta.*const3)*ve_wo+mgOK*ac_wo+((2*neta/delt).*const2)*de_vo+neta.*const2*ve_vo-fvg;  
   
   for linha=1:Nnos
       fcell{linha,1}=[f1(linha,1);f2(linha,1)];
   end
   vetor_F=cell2mat(fcell);
   vetor_F=vetor_F+(1*aci(:,inst)); %OBS: na validação entrávamos direto com F, então "vetor_F=vetor_F+(1*aci(:,inst))", aci era a força já (34)
   
   restart=4; tol=1e-8; maxit=400;
   [resul,flag,relres,iter,resvec] = gmres(matriz_A(ide_L,ide_L),vetor_F(ide_L),restart,tol,maxit,L,U);
   converg(inst,1)=flag; 
   flag  % flag=0: houve convergência para a tol desejada! flag=1: não houve convergência com o maxit colocado 
   
   deslocamento(ide_L,inst)=resul;
   ve_x(ide_Lno,inst)=(2/delt).*(resul(1:2:length(resul))-(de_vo(ide_Lno,1)+delt*0.5*ve_vo(ide_Lno,1)));
   ve_y(ide_Lno,inst)=(2/delt).*(resul(2:2:length(resul))-(de_wo(ide_Lno,1)+delt*0.5*ve_wo(ide_Lno,1)));   
   ac_x(ide_Lno,inst)=(4/delt^2)*(resul(1:2:length(resul))-(de_vo(ide_Lno,1)+delt*ve_vo(ide_Lno,1)+0.25*(delt^2)*ac_vo(ide_Lno,1)));
   ac_y(ide_Lno,inst)=(4/delt^2)*(resul(2:2:length(resul))-(de_wo(ide_Lno,1)+delt*ve_wo(ide_Lno,1)+0.25*(delt^2)*ac_wo(ide_Lno,1))); 
  
   de_vo(ide_Lno,1)=resul(1:2:length(resul)); de_wo(ide_Lno,1)=resul(2:2:length(resul)); 
   ve_vo=ve_x(:,inst); ve_wo=ve_y(:,inst); ac_vo=ac_x(:,inst); ac_wo=ac_y(:,inst);
   contt=contt+1
   
   
   % Condicional para verificar se já chegou no regime permanente
   if inst >= 500 && mod(inst,500) == 0
       disp('Testando o regime')
       [regime, idc] = regimeCheck(ac_y(floor(Gsup/2),ini:contt-1),Fs,we,5,0.1);
       if strcmp(regime,'steady')
           disp(regime)
           break
       end
       ini = contt;
   end

end

% Construindo os vetores de velocidade e aceleração
% for coluna=1:contt-1
%    for linha=1:Nnos
%       velocidade{linha,coluna}=[ve_x(linha,coluna);ve_y(linha,coluna)];
%       aceleracao{linha,coluna}=[ac_x(linha,coluna);ac_y(linha,coluna)];
%    end
% end
% velocidade=cell2mat(velocidade);
% aceleracao=cell2mat(aceleracao);

% Forma otimizada
velocidade = zeros(Nnos,(contt-1)*2);
velocidade(1:2:end) = ve_x(:,1:contt-1);
velocidade(2:2:end) = ve_y(:,1:contt-1);

aceleracao = zeros(Nnos,(contt-1)*2);
aceleracao(1:2:end) = ac_x(:,1:contt-1);
aceleracao(2:2:end) = ac_y(:,1:contt-1);


%Gsup=local;
% figure(3)
% plot(t,deslocamento(Gsup,:))
% xlim([0.01 tf])
% title('Deslocamento')
% xlabel('Tempo [s]','FontSize',10)
% ylabel('Deslocamento [m]','FontSize',10)
% 
% figure(4)
% plot(t,velocidade(Gsup,:))
% xlim([0.02 tf])
% title('Velocidade')
% xlabel('Tempo [s]','FontSize',10)
% ylabel('Velocidade [m/s]','FontSize',10)
% 
% figure(5)
% plot(t,aceleracao(Gsup,:))
% xlim([0.04 tf])
% title('Aceleração')
% xlabel('Tempo [s]','FontSize',10)
% ylabel('Aceleração [m/s^2]','FontSize',10)
% 
% figure(6)
% plot(t,aceleracao(local,:),'r',t,aceleracao(Gsup,:),'b')
% grid on
% xlim([0.8 0.86])
% title('Comparação entre a aceleração do nó excitado da base e o nó sup central')
% xlabel('Tempo [s]','FontSize',10)
% ylabel('Aceleração [m/s^2]','FontSize',10)

% Checando se a aceleração do nó da fonte está ok (esse inicio de 3000 é para tirar o intervalo de acomodação)
% RMS_F=rms(aceleracao(local,1:contt:length(aceleracao(local,:)))) 
% RMS_S1=rms(aceleracao(Gsup,3000:length(aceleracao(Gsup,:))))
% Atenuacao=((RMS_F-RMS_S1)/RMS_F) % Fonte x Superfície

% NFFT_ac = 2^nextpow2(length(aceleracao(local,2000:length(aceleracao(local,:)))));
% FFT_S1 = fft(aceleracao(Gsup,2000:length(aceleracao(local,:))),NFFT_ac)/length(aceleracao(Gsup,2000:length(aceleracao(local,:))));
% FFT_F = fft(aceleracao(local,2000:length(aceleracao(local,:))),NFFT_ac)/length(aceleracao(local,2000:length(aceleracao(local,:))));
% ff=Fs/2*linspace(0,1,NFFT_ac/2+1);
% 
% Mmo = 0;
% Mmax = 0;
% 
% for i=5:length(ff)
%  FFT_sinal1_A_ac=abs(FFT_S1(i -1));
%  FFT_sinal1_B_ac=abs(FFT_S1(i));
%  if (FFT_sinal1_B_ac > FFT_sinal1_A_ac)
%    Mmo = FFT_sinal1_B_ac ;
%    if (Mmo>= Mmax)
%      Mmax = Mmo ;
%      fase_sinal1_ac = phase(FFT_S1(i));
%    end
%  end
% end
% 
% Mmo = 0;
% Mmax = 0;
% 
% for i = 5:length(ff)
%  FFT_sinal2_A_ac = abs(FFT_F(i-1));
%  FFT_sinal2_B_ac = abs(FFT_F(i));
%  if (FFT_sinal2_B_ac > FFT_sinal2_A_ac)
%    Mmo = FFT_sinal2_B_ac;
%    if (Mmo>= Mmax)
%      Mmax = Mmo ;
%      fase_sinal2_ac = phase (FFT_F(i));
%    end
%  end
% end
% 
% defasagem_ac=abs(fase_sinal1_ac-fase_sinal2_ac)*(180/pi)
% -------------------------------------------------------------------------