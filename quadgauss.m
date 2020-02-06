%--------------------------------------------------------------------------   
%                  IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%               Universidade Estadual Paulista, FEIS/UNESP

%                              quadgauss
%      Implementação do Método numérico para resolução das integrais

% ------------------------ QUADRATURA DE GAUSS ----------------------------

function [P,W]=quadgauss(np)

if np==1,
P(1)=0;
W(1)=2;
end

if np==2,
P(1)=-1/sqrt(3);
P(2)=+1/sqrt(3);
W(1)=1;
W(2)=1;
end

if np==3,
P(1)=-sqrt(3)/sqrt(5);
P(2)=0;
P(3)=sqrt(3)/sqrt(5);
W(1)=5/9;
W(2)=8/9;
W(3)=5/9;
end

if np==4,
P(1)=-0.8611363116;
P(2)=-0.3399810436;
P(3)=0.3399810436;
P(4)=0.8611363116;
W(1)=0.3478548451;
W(2)=0.6521451549;
W(3)=0.6521451549;
W(4)=0.3478548451;
end