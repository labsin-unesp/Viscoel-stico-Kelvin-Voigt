%--------------------------------------------------------------------------   
%                  IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%               Universidade Estadual Paulista, FEIS/UNESP

%                               Jacobiano
%                 Matriz de transformação das variáveis

%--------------------------- MATRIZ JACOBIANA -----------------------------

function [J11i,J12i,J21i,J22i,detJ]=jacobiano(x11,x21,x31,x41,x12,x22,x32,x42,P)

np=length(P);
for n1=1:np,
    for n2=1:np,
        s1=P(n1);
        s2=P(n2);
        J11(n1,n2)=(1/4)*((-1+s2)*x11+(+1-s2)*x21+(+1+s2)*x31+(-1-s2)*x41);
        J12(n1,n2)=(1/4)*((-1+s1)*x11+(-1-s1)*x21+(+1+s1)*x31+(+1-s1)*x41);
        J21(n1,n2)=(1/4)*((-1+s2)*x12+(+1-s2)*x22+(+1+s2)*x32+(-1-s2)*x42);
        J22(n1,n2)=(1/4)*((-1+s1)*x12+(-1-s1)*x22+(+1+s1)*x32+(+1-s1)*x42);
    end
end

detJ=J11.*J22-J12.*J21;
J11i=J22./detJ;
J12i=-J12./detJ;
J21i=-J21./detJ;
J22i=J11./detJ;
