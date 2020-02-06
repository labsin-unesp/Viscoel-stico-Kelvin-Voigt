%--------------------------------------------------------------------------   
%                 IMPLEMENTA��O NUM�RICA - DISSERTA��O
%               Universidade Estadual Paulista, FEIS/UNESP

%                               restri��o:
%      Constroi o vetor com os graus de liberdade restritos do problema
% 
%  No caso de um duplo engaste, o vetor ser� formado por todos os graus de
%             liberdade dos n�s das extremidades laterais.

%OBS: Quanto a numera��o - O n� 1, por exemplo, tem os GL 1 e 2. 1 refere a x e 2 a y. 
%                          O n� 20 tem os GL 39 e 40, 39 refere a x e 40 a y.

%-------------------------- VETOR DOS GL LIVRES ---------------------------

% ENGASTADO DOS LADOS:
function [ide_P]=condcont(Ndivx,Ndivy)

Nx=Ndivx+1;                % N�mero de n�s em x
Ny=Ndivy+1;                % N�mero de n�s em y
ide_P=[];                  % Armazena os n�s livres
cc=1;
cl=0;

while cc<=Ny
    ide_P=[ide_P,1+cl,2+cl,(2*Nx-1)+cl,(2*Nx)+cl];
    cc=cc+1;
    cl=cl+2*Nx;
end    
end