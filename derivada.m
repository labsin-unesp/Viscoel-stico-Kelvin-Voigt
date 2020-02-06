% Esta function calcula a derivada de uma fun��o dada por um conjunto de pontos discretos.

function dx=derivada(x,y)
n=length(x);
dx(1)=(y(2)-y(1))/(x(2)-x(1));     % diferen�as finitas progressiva
for i=2:n-1
    dx(i)=(y(i+1)-y(i-1))/(x(i+1)-x(i-1));
end
dx(n)=(y(n)-y(n-1))/(x(n)-x(n-1)); % diferen�as finitas regressiva
