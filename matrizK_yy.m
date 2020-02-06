%--------------------------------------------------------------------------   
%                 IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%             Universidade Estadual Paulista, FEIS/UNESP

%                  - Calculo da "matriz de rigidez" - 
%    * Parcela referente a derivada da funções de forma em relação a y *
%(Lembrando que estamos em um plano x-y - na formulação teórica fizemos y-z)
% x=y; y=z

%-------------------------- MATRIZ DE RIGIDEZ -----------------------------

function k_yy=matrizK_yy(h,detJ,B21,B22,B23,B24,B25,B26,B27,B28,np)

[P,W]=quadgauss(np);
k11=0;
k12=0;
k13=0;
k14=0;
k15=0;
k16=0;
k17=0;
k18=0;
k21=0;
k22=0;
k23=0;
k24=0;
k25=0;
k26=0;
k27=0;
k28=0;
k31=0;
k32=0;
k33=0;
k34=0;
k35=0;
k36=0;
k37=0;
k38=0;
k41=0;
k42=0;
k43=0;
k44=0;
k45=0;
k46=0;
k47=0;
k48=0;
k51=0;
k52=0;
k53=0;
k54=0;
k55=0;
k56=0;
k57=0;
k58=0;
k61=0;
k62=0;
k63=0;
k64=0;
k65=0;
k66=0;
k67=0;
k68=0;
k71=0;
k72=0;
k73=0;
k74=0;
k75=0;
k76=0;
k77=0;
k78=0;
k81=0;
k82=0;
k83=0;
k84=0;
k85=0;
k86=0;
k87=0;
k88=0;

for n1=1:np,
    for n2=1:np,
        k11=k11+W(n1)*W(n2)*(B21(n1,n2)*B21(n1,n2))*detJ(n1,n2);
        k12=k12+W(n1)*W(n2)*(B22(n1,n2)*B21(n1,n2))*detJ(n1,n2);
        k13=k13+W(n1)*W(n2)*(B23(n1,n2)*B21(n1,n2))*detJ(n1,n2);
        k14=k14+W(n1)*W(n2)*(B24(n1,n2)*B21(n1,n2))*detJ(n1,n2);
        k15=k15+W(n1)*W(n2)*(B25(n1,n2)*B21(n1,n2))*detJ(n1,n2);
        k16=k16+W(n1)*W(n2)*(B26(n1,n2)*B21(n1,n2))*detJ(n1,n2);
        k17=k17+W(n1)*W(n2)*(B27(n1,n2)*B21(n1,n2))*detJ(n1,n2);
        k18=k18+W(n1)*W(n2)*(B28(n1,n2)*B21(n1,n2))*detJ(n1,n2);
        
        k21=k21+W(n1)*W(n2)*(B21(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        k22=k22+W(n1)*W(n2)*(B22(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        k23=k23+W(n1)*W(n2)*(B23(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        k24=k24+W(n1)*W(n2)*(B24(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        k25=k25+W(n1)*W(n2)*(B25(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        k26=k26+W(n1)*W(n2)*(B26(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        k27=k27+W(n1)*W(n2)*(B27(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        k28=k28+W(n1)*W(n2)*(B28(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        
        k31=k31+W(n1)*W(n2)*(B21(n1,n2)*B23(n1,n2))*detJ(n1,n2);
        k32=k32+W(n1)*W(n2)*(B22(n1,n2)*B23(n1,n2))*detJ(n1,n2);
        k33=k33+W(n1)*W(n2)*(B23(n1,n2)*B23(n1,n2))*detJ(n1,n2);
        k34=k34+W(n1)*W(n2)*(B24(n1,n2)*B23(n1,n2))*detJ(n1,n2);
        k35=k35+W(n1)*W(n2)*(B25(n1,n2)*B23(n1,n2))*detJ(n1,n2);
        k36=k36+W(n1)*W(n2)*(B26(n1,n2)*B23(n1,n2))*detJ(n1,n2);
        k37=k37+W(n1)*W(n2)*(B27(n1,n2)*B23(n1,n2))*detJ(n1,n2);
        k38=k38+W(n1)*W(n2)*(B28(n1,n2)*B23(n1,n2))*detJ(n1,n2);
        
        k41=k41+W(n1)*W(n2)*(B21(n1,n2)*B24(n1,n2))*detJ(n1,n2);
        k42=k42+W(n1)*W(n2)*(B22(n1,n2)*B24(n1,n2))*detJ(n1,n2);
        k43=k43+W(n1)*W(n2)*(B23(n1,n2)*B24(n1,n2))*detJ(n1,n2);
        k44=k44+W(n1)*W(n2)*(B24(n1,n2)*B24(n1,n2))*detJ(n1,n2);
        k45=k45+W(n1)*W(n2)*(B25(n1,n2)*B24(n1,n2))*detJ(n1,n2);
        k46=k46+W(n1)*W(n2)*(B26(n1,n2)*B24(n1,n2))*detJ(n1,n2);
        k47=k47+W(n1)*W(n2)*(B27(n1,n2)*B24(n1,n2))*detJ(n1,n2);
        k48=k48+W(n1)*W(n2)*(B28(n1,n2)*B24(n1,n2))*detJ(n1,n2);
        
        k51=k51+W(n1)*W(n2)*(B21(n1,n2)*B25(n1,n2))*detJ(n1,n2);
        k52=k52+W(n1)*W(n2)*(B22(n1,n2)*B25(n1,n2))*detJ(n1,n2);
        k53=k53+W(n1)*W(n2)*(B23(n1,n2)*B25(n1,n2))*detJ(n1,n2);
        k54=k54+W(n1)*W(n2)*(B24(n1,n2)*B25(n1,n2))*detJ(n1,n2);
        k55=k55+W(n1)*W(n2)*(B25(n1,n2)*B25(n1,n2))*detJ(n1,n2);
        k56=k56+W(n1)*W(n2)*(B26(n1,n2)*B25(n1,n2))*detJ(n1,n2);
        k57=k57+W(n1)*W(n2)*(B27(n1,n2)*B25(n1,n2))*detJ(n1,n2);
        k58=k58+W(n1)*W(n2)*(B28(n1,n2)*B25(n1,n2))*detJ(n1,n2);
        
        k61=k61+W(n1)*W(n2)*(B21(n1,n2)*B26(n1,n2))*detJ(n1,n2);
        k62=k62+W(n1)*W(n2)*(B22(n1,n2)*B26(n1,n2))*detJ(n1,n2);
        k63=k63+W(n1)*W(n2)*(B23(n1,n2)*B26(n1,n2))*detJ(n1,n2);
        k64=k64+W(n1)*W(n2)*(B24(n1,n2)*B26(n1,n2))*detJ(n1,n2);
        k65=k65+W(n1)*W(n2)*(B25(n1,n2)*B26(n1,n2))*detJ(n1,n2);
        k66=k66+W(n1)*W(n2)*(B26(n1,n2)*B26(n1,n2))*detJ(n1,n2);
        k67=k67+W(n1)*W(n2)*(B27(n1,n2)*B26(n1,n2))*detJ(n1,n2);
        k68=k68+W(n1)*W(n2)*(B28(n1,n2)*B26(n1,n2))*detJ(n1,n2);
        
        k71=k71+W(n1)*W(n2)*(B21(n1,n2)*B27(n1,n2))*detJ(n1,n2);
        k72=k72+W(n1)*W(n2)*(B22(n1,n2)*B27(n1,n2))*detJ(n1,n2);
        k73=k73+W(n1)*W(n2)*(B23(n1,n2)*B27(n1,n2))*detJ(n1,n2);
        k74=k74+W(n1)*W(n2)*(B24(n1,n2)*B27(n1,n2))*detJ(n1,n2);
        k75=k75+W(n1)*W(n2)*(B25(n1,n2)*B27(n1,n2))*detJ(n1,n2);
        k76=k76+W(n1)*W(n2)*(B26(n1,n2)*B27(n1,n2))*detJ(n1,n2);
        k77=k77+W(n1)*W(n2)*(B27(n1,n2)*B27(n1,n2))*detJ(n1,n2);
        k78=k78+W(n1)*W(n2)*(B28(n1,n2)*B27(n1,n2))*detJ(n1,n2);
        
        k81=k81+W(n1)*W(n2)*(B21(n1,n2)*B28(n1,n2))*detJ(n1,n2);
        k82=k82+W(n1)*W(n2)*(B22(n1,n2)*B28(n1,n2))*detJ(n1,n2);
        k83=k83+W(n1)*W(n2)*(B23(n1,n2)*B28(n1,n2))*detJ(n1,n2);
        k84=k84+W(n1)*W(n2)*(B24(n1,n2)*B28(n1,n2))*detJ(n1,n2);
        k85=k85+W(n1)*W(n2)*(B25(n1,n2)*B28(n1,n2))*detJ(n1,n2);
        k86=k86+W(n1)*W(n2)*(B26(n1,n2)*B28(n1,n2))*detJ(n1,n2);
        k87=k87+W(n1)*W(n2)*(B27(n1,n2)*B28(n1,n2))*detJ(n1,n2);
        k88=k88+W(n1)*W(n2)*(B28(n1,n2)*B28(n1,n2))*detJ(n1,n2);   
    end 
end

k_yy=h*[
k11,k12,k13,k14,k15,k16,k17,k18
k21,k22,k23,k24,k25,k26,k27,k28
k31,k32,k33,k34,k35,k36,k37,k38
k41,k42,k43,k44,k45,k46,k47,k48
k51,k52,k53,k54,k55,k56,k57,k58
k61,k62,k63,k64,k65,k66,k67,k68
k71,k72,k73,k74,k75,k76,k77,k78
k81,k82,k83,k84,k85,k86,k87,k88];
