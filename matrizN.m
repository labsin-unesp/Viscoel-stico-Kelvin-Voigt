%--------------------------------------------------------------------------   
%                  IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%               Universidade Estadual Paulista, FEIS/UNESP

%                                matrizN
%                   Matriz das Funções interpoladoras

%------------------------------ MATRIZ N ----------------------------------


function [N11,N12,N13,N14,N15,N16,N17,N18,N21,N22,N23,N24,N25,N26,N27,N28]=matrizN(np)

[P,W]=quadgauss(np);
for n1=1:np,
    for n2=1:np,
        s1=P(n1);
        s2=P(n2);
        N11(n1,n2)=(1-s1)*(1-s2)/4;
        N12(n1,n2)=0;
        N13(n1,n2)=(1+s1)*(1-s2)/4;
        N14(n1,n2)=0;
        N15(n1,n2)=(1+s1)*(1+s2)/4;
        N16(n1,n2)=0;
        N17(n1,n2)=(1-s1)*(1+s2)/4;
        N18(n1,n2)=0;
        N21(n1,n2)=0;
        N22(n1,n2)=(1-s1)*(1-s2)/4;
        N23(n1,n2)=0;
        N24(n1,n2)=(1+s1)*(1-s2)/4;
        N25(n1,n2)=0;
        N26(n1,n2)=(1+s1)*(1+s2)/4;
        N27(n1,n2)=0;
        N28(n1,n2)=(1-s1)*(1+s2)/4;
    end
end
