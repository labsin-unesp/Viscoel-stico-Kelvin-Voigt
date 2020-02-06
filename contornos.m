%--------------------------------------------------------------------------   
%                  IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%               Universidade Estadual Paulista, FEIS/UNESP

%                  - Calculo nos contornos gama1 e gama2 - 
%    
%(Lembrando que estamos em um plano x-y - na formulação teórica fizemos y-z)
% x=y; y=z

%----------------------------- CONTORNOS ----------------------------------

function gama=contornos(h,detJ,N11,N12,N13,N14,N15,N16,N17,N18,B21,B22,B23,B24,B25,B26,B27,B28,np)

[P,W]=quadgauss(np);
C11=0;
C12=0;
C13=0;
C14=0;
C15=0;
C16=0;
C17=0;
C18=0;
C21=0;
C22=0;
C23=0;
C24=0;
C25=0;
C26=0;
C27=0;
C28=0;
C31=0;
C32=0;
C33=0;
C34=0;
C35=0;
C36=0;
C37=0;
C38=0;
C41=0;
C42=0;
C43=0;
C44=0;
C45=0;
C46=0;
C47=0;
C48=0;
C51=0;
C52=0;
C53=0;
C54=0;
C55=0;
C56=0;
C57=0;
C58=0;
C61=0;
C62=0;
C63=0;
C64=0;
C65=0;
C66=0;
C67=0;
C68=0;
C71=0;
C72=0;
C73=0;
C74=0;
C75=0;
C76=0;
C77=0;
C78=0;
C81=0;
C82=0;
C83=0;
C84=0;
C85=0;
C86=0;
C87=0;
C88=0;

for n1=1:np,
    for n2=1:np,
        C11=C11+W(n1)*W(n2)*(N11(n1,n2)*B21(n1,n2))*detJ(n1,n2); 
        C12=C12+W(n1)*W(n2)*(N11(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        C13=C13+W(n1)*W(n2)*(N11(n1,n2)*B23(n1,n2))*detJ(n1,n2); 
        C14=C14+W(n1)*W(n2)*(N11(n1,n2)*B24(n1,n2))*detJ(n1,n2); 
        C15=C15+W(n1)*W(n2)*(N11(n1,n2)*B25(n1,n2))*detJ(n1,n2); 
        C16=C16+W(n1)*W(n2)*(N11(n1,n2)*B26(n1,n2))*detJ(n1,n2); 
        C17=C17+W(n1)*W(n2)*(N11(n1,n2)*B27(n1,n2))*detJ(n1,n2); 
        C18=C18+W(n1)*W(n2)*(N11(n1,n2)*B28(n1,n2))*detJ(n1,n2); 
        
        C21=C21+W(n1)*W(n2)*(N12(n1,n2)*B21(n1,n2))*detJ(n1,n2); 
        C22=C22+W(n1)*W(n2)*(N12(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        C23=C23+W(n1)*W(n2)*(N12(n1,n2)*B23(n1,n2))*detJ(n1,n2); 
        C24=C24+W(n1)*W(n2)*(N12(n1,n2)*B24(n1,n2))*detJ(n1,n2); 
        C25=C25+W(n1)*W(n2)*(N12(n1,n2)*B25(n1,n2))*detJ(n1,n2); 
        C26=C26+W(n1)*W(n2)*(N12(n1,n2)*B26(n1,n2))*detJ(n1,n2); 
        C27=C27+W(n1)*W(n2)*(N12(n1,n2)*B27(n1,n2))*detJ(n1,n2); 
        C28=C28+W(n1)*W(n2)*(N12(n1,n2)*B28(n1,n2))*detJ(n1,n2); 
        
        C31=C31+W(n1)*W(n2)*(N13(n1,n2)*B21(n1,n2))*detJ(n1,n2); 
        C32=C32+W(n1)*W(n2)*(N13(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        C33=C33+W(n1)*W(n2)*(N13(n1,n2)*B23(n1,n2))*detJ(n1,n2); 
        C34=C34+W(n1)*W(n2)*(N13(n1,n2)*B24(n1,n2))*detJ(n1,n2); 
        C35=C35+W(n1)*W(n2)*(N13(n1,n2)*B25(n1,n2))*detJ(n1,n2); 
        C36=C36+W(n1)*W(n2)*(N13(n1,n2)*B26(n1,n2))*detJ(n1,n2); 
        C37=C37+W(n1)*W(n2)*(N13(n1,n2)*B27(n1,n2))*detJ(n1,n2); 
        C38=C38+W(n1)*W(n2)*(N13(n1,n2)*B28(n1,n2))*detJ(n1,n2); 
        
        C41=C41+W(n1)*W(n2)*(N14(n1,n2)*B21(n1,n2))*detJ(n1,n2); 
        C42=C42+W(n1)*W(n2)*(N14(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        C43=C43+W(n1)*W(n2)*(N14(n1,n2)*B23(n1,n2))*detJ(n1,n2); 
        C44=C44+W(n1)*W(n2)*(N14(n1,n2)*B24(n1,n2))*detJ(n1,n2); 
        C45=C45+W(n1)*W(n2)*(N14(n1,n2)*B25(n1,n2))*detJ(n1,n2); 
        C46=C46+W(n1)*W(n2)*(N14(n1,n2)*B26(n1,n2))*detJ(n1,n2); 
        C47=C47+W(n1)*W(n2)*(N14(n1,n2)*B27(n1,n2))*detJ(n1,n2); 
        C48=C48+W(n1)*W(n2)*(N14(n1,n2)*B28(n1,n2))*detJ(n1,n2); 
        
        C51=C51+W(n1)*W(n2)*(N15(n1,n2)*B21(n1,n2))*detJ(n1,n2); 
        C52=C52+W(n1)*W(n2)*(N15(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        C53=C53+W(n1)*W(n2)*(N15(n1,n2)*B23(n1,n2))*detJ(n1,n2); 
        C54=C54+W(n1)*W(n2)*(N15(n1,n2)*B24(n1,n2))*detJ(n1,n2); 
        C55=C55+W(n1)*W(n2)*(N15(n1,n2)*B25(n1,n2))*detJ(n1,n2); 
        C56=C56+W(n1)*W(n2)*(N15(n1,n2)*B26(n1,n2))*detJ(n1,n2); 
        C57=C57+W(n1)*W(n2)*(N15(n1,n2)*B27(n1,n2))*detJ(n1,n2); 
        C58=C58+W(n1)*W(n2)*(N15(n1,n2)*B28(n1,n2))*detJ(n1,n2); 
        
        C61=C61+W(n1)*W(n2)*(N16(n1,n2)*B21(n1,n2))*detJ(n1,n2); 
        C62=C62+W(n1)*W(n2)*(N16(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        C63=C63+W(n1)*W(n2)*(N16(n1,n2)*B23(n1,n2))*detJ(n1,n2); 
        C64=C64+W(n1)*W(n2)*(N16(n1,n2)*B24(n1,n2))*detJ(n1,n2); 
        C65=C65+W(n1)*W(n2)*(N16(n1,n2)*B25(n1,n2))*detJ(n1,n2); 
        C66=C66+W(n1)*W(n2)*(N16(n1,n2)*B26(n1,n2))*detJ(n1,n2); 
        C67=C67+W(n1)*W(n2)*(N16(n1,n2)*B27(n1,n2))*detJ(n1,n2); 
        C68=C68+W(n1)*W(n2)*(N16(n1,n2)*B28(n1,n2))*detJ(n1,n2); 
        
        C71=C71+W(n1)*W(n2)*(N17(n1,n2)*B21(n1,n2))*detJ(n1,n2); 
        C72=C72+W(n1)*W(n2)*(N17(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        C73=C73+W(n1)*W(n2)*(N17(n1,n2)*B23(n1,n2))*detJ(n1,n2); 
        C74=C74+W(n1)*W(n2)*(N17(n1,n2)*B24(n1,n2))*detJ(n1,n2); 
        C75=C75+W(n1)*W(n2)*(N17(n1,n2)*B25(n1,n2))*detJ(n1,n2); 
        C76=C76+W(n1)*W(n2)*(N17(n1,n2)*B26(n1,n2))*detJ(n1,n2); 
        C77=C77+W(n1)*W(n2)*(N17(n1,n2)*B27(n1,n2))*detJ(n1,n2); 
        C78=C78+W(n1)*W(n2)*(N17(n1,n2)*B28(n1,n2))*detJ(n1,n2); 
        
        C81=C81+W(n1)*W(n2)*(N18(n1,n2)*B21(n1,n2))*detJ(n1,n2); 
        C82=C82+W(n1)*W(n2)*(N18(n1,n2)*B22(n1,n2))*detJ(n1,n2);
        C83=C83+W(n1)*W(n2)*(N18(n1,n2)*B23(n1,n2))*detJ(n1,n2); 
        C84=C84+W(n1)*W(n2)*(N18(n1,n2)*B24(n1,n2))*detJ(n1,n2); 
        C85=C85+W(n1)*W(n2)*(N18(n1,n2)*B25(n1,n2))*detJ(n1,n2); 
        C86=C86+W(n1)*W(n2)*(N18(n1,n2)*B26(n1,n2))*detJ(n1,n2); 
        C87=C87+W(n1)*W(n2)*(N18(n1,n2)*B27(n1,n2))*detJ(n1,n2); 
        C88=C88+W(n1)*W(n2)*(N18(n1,n2)*B28(n1,n2))*detJ(n1,n2); 
    end  
end

gama=h*[
C11,C12,C13,C14,C15,C16,C17,C18
C21,C22,C23,C24,C25,C26,C27,C28
C31,C32,C33,C34,C35,C36,C37,C38
C41,C42,C43,C44,C45,C46,C47,C48
C51,C52,C53,C54,C55,C56,C57,C58
C61,C62,C63,C64,C65,C66,C67,C68
C71,C72,C73,C74,C75,C76,C77,C78
C81,C82,C83,C84,C85,C86,C87,C88];