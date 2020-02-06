%--------------------------------------------------------------------------   
%                 IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%              Universidade Estadual Paulista, FEIS/UNESP

%                                ice_G 
%           Matriz com os graus de liberdade de cada elemento
%                          Obs: 2gld por nó

%--------------------------------------------------------------------------

function [ice_G]=iceG(Nele,ice_N)

for e=1:1:Nele,
    ice_G(e,1)=ice_N(e,1)*2-1;
    ice_G(e,2)=ice_N(e,1)*2;
    ice_G(e,3)=ice_N(e,2)*2-1;
    ice_G(e,4)=ice_N(e,2)*2;
    ice_G(e,5)=ice_N(e,3)*2-1;
    ice_G(e,6)=ice_N(e,3)*2;
    ice_G(e,7)=ice_N(e,4)*2-1;
    ice_G(e,8)=ice_N(e,4)*2;
end
