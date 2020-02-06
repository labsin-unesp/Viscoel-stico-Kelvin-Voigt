%--------------------------------------------------------------------------   
%                  IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%               Universidade Estadual Paulista, FEIS/UNESP

%     Implementando a parcela referente á força Gravitacional atuante 

%-------------------------- FORÇA GRAVITACIONAL ---------------------------


function fv=gravidade(N11,N13,N15,N17,p,h,detJ,np)
%b=[0;-1]; % componentes das forças de volume
[P,W]=quadgauss(np);
fv11=0;
fv12=0;
fv21=0;
fv22=0;
fv31=0;
fv32=0;
fv41=0;
fv42=0;
for n1=1:np
    for n2=1:np
        fv11=fv11+0;
        fv12=fv12+W(n1)*W(n2)*detJ(n1,n2)*N11(n1,n2);
        fv21=fv21+0;
        fv22=fv22+W(n1)*W(n2)*detJ(n1,n2)*N13(n1,n2);
        fv31=fv31+0;
        fv32=fv32+W(n1)*W(n2)*detJ(n1,n2)*N15(n1,n2);
        fv41=fv41+0;
        fv42=fv42+W(n1)*W(n2)*detJ(n1,n2)*N17(n1,n2);
    end
end
fv=p(1)*h(1)*(-9.81)*[fv12;fv22;fv32;fv42];
