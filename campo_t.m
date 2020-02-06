function Z = campo_t(Ndivx,Ndivy,aceleracao,t)
k = 1;

for i = 1:Ndivx

    for j = 1:Ndivy
        
        disp([i,j,k])
        
        Zx(i,j) = aceleracao(k,t);
        Zy(i,j) = aceleracao(k+1,t);
        
        k = k + 2;
        
        Z(i,j) = sqrt(Zx(i,j)^2 +Zy(i,j)^2);
    end
end

end