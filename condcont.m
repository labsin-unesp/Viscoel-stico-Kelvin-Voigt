%--------------------------------------------------------------------------   
%                 IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%               Universidade Estadual Paulista, FEIS/UNESP

%                               condcont:
%      Constroi o vetor com os graus de liberdade livres do problema
% 
% No caso abordado assumiremos um duplo engaste. as extremidades laterais
%   estarão com os graus de liberdade restritos tanto em x quanto em y

%OBS: Quanto a numeração - O nó 1, por exemplo, tem os GL 1 e 2. 1 refere a x e 2 a y. 
%                          O nó 20 tem os GL 39 e 40, 39 refere a x e 40 a y.

%-------------------------- VETOR DOS GL LIVRES ---------------------------

% ENGASTADO DOS LADOS:
function [ide_L]=condcont(Ndivx,Ndivy)

Nx=Ndivx+1;                % Número de nós em x
Ny=Ndivy+1;                % Número de nós em y
ide_L=[];                  % Armazena os nós livres
cc=1;
cl=0;

while cc<=Ny
    ide_L=[ide_L,(2*(2+cl)-1):(2*((Nx-1)+cl))];
    cc=cc+1;
    cl=cl+Nx;
end    
end

% -------------------------------------------------------------------------
% ENGASTADO NA BASE E PRESO EM X DOS LADOS:

% function [ide_L]=condcont(Ndivx,Ndivy)
% ide_L=[];
% local=2*(floor((Ndivx+2)/2));            % liberando o GLy do nó no qual encontra-se a excitação
% ide_L=[ide_L,local-10,local-8,local-6,local-4,local-2,local,local+2,local+4,local+6,local+8,local+10];
% for j=1:Ndivy
% ide_L=[ide_L,2*(((Ndivx+1)*j)+1):2*((Ndivx+1)*(j+1))-2,2*(Ndivx+1)*(j+1)];
% end
