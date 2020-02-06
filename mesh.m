%--------------------------------------------------------------------------   
%                  IMPLEMENTAÇÃO NUMÉRICA - DISSERTAÇÃO
%               Universidade Estadual Paulista, FEIS/UNESP

%                                  mesh
%             Construção automática da malha a ser utilizada
%                Obs: Elementos quadrilaterais de 4nós 

%-------------------------- GERADOR DE MALHA ------------------------------

function [x,y,ice_N]=mesh(xv,yv,Ndivx,Ndivy,Nele)
h=0;
a=1;
k=1;
while a<=Ndivy+1
    while k<=(a*(Ndivx+1))
    x(k)= 0 + h*(xv/(Ndivx));
    h=h+1;
    k=k+1;
    end
h=0;
a=a+1;
end
j=1;
b=1;
while b<=Ndivy+1
    while j<=(b*(Ndivx+1))
    y(j)= 0 + h*(yv/(Ndivy));
    j=j+1;
    end
h=h+1;
b=b+1;
end
d=1;
r=1;
q=1;
d=1;
i=0;
while d<=Nele
    while r<=(q*Ndivx+i)
        ice_N(d,1)=r;
        ice_N(d,2)=ice_N(d,1)+1;
        ice_N(d,3)=ice_N(d,2)+(Ndivx+1);
        ice_N(d,4)=ice_N(d,3)-1;
        d=d+1;
        r=r+1;
    end
    i=i+1;
    r=r+1;
    q=q+1;
end
end

