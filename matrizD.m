%--------------------------------------------------------------------------   
%                 IMPLEMENTA��O NUM�RICA - DISSERTA��O
%              Universidade Estadual Paulista, FEIS/UNESP

%                   [D] Matriz constitutiva el�stica
%                  {Tens�es}=[D]*{def + neta*taxadef}
%--------------------------------------------------------------------------

% PARA O ESTADO PLANO DE TENS�ES:
function D=matrizD(E,v)

G=E/(2*(1+v));  % M�dulo cisalhante

D(1,1) = E/(1-v^2);
D(1,2) = (E*v)/(1-v^2);
D(1,3) = 0;
D(2,1) = (E*v)/(1-v^2);
D(2,2) = E/(1-v^2);
D(2,3) = 0;
D(3,1) = 0;
D(3,2) = 0;
D(3,3) = G;

% OBS. Def_z= -(v/E)*(Ten_x + Ten_y) = -((v(1+v))/(1-v^2))*(Def_x +Def_y);  (terceiro eixo)


% -------------------------------------------------------------------------
% PARA O ESTADO PLANO DE DEFORMA��O:

% function D=matrizD(E,v)
% 
% bulk=E/(3*(1-2*v));           % K, M�dulo de expans�o volum�trica (Bulk modulus)
% G=E/(2*(1+v));                % M�dulo cisalhante
% 
% D(1,1) = bulk+(4*G/3);
% D(1,2) = bulk-(2*G/3);
% D(1,3) = 0;
% D(2,1) = bulk-(2*G/3);
% D(2,2) = bulk+(4*G/3);
% D(2,3) = 0;
% D(3,1) = 0;
% D(3,2) = 0;
% D(3,3) = G;

% OBS. Ten_z= v*(Ten_x + Ten_y);
