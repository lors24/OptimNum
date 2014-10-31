function [ x0, fx, iter ] = MiSimplex(c, A, b, x0 )

%Haces todo columnas
c = c(:);
b = b(:);
x0 = x0(:);

[m,n] = size(A);

sol = 0;
iter = 0;

indb = find(x0);                                                            %encuentra las variables b�sicas del v�rtice (�ndices) // vector columna
indnb = find(x0 == 0);                                                      %encuentra las variables no b�sicas del v�rtice (�ndices)  // vector columna


while (sol == 0)  %&& (A*(x0) == b)
iter = iter+1;   

%si el v�rtice es no degenerado
%completa los indices de variables b�sicas
   
if length(indb) < m
i = 1;  
while ( length(indb) < m  )
    aux = [indb;indnb(i)];
    AUX = A(:,aux);
    if  rank(AUX) == length(aux)
        indb = aux;                                                         %nuevo vector �ndices de variables b�sicas // vector columna
        indnb(i) = 0;
    end
    i = i+1;
end

aux2 = find(indnb);
indnb = indnb(aux2);
indnb = indnb(:);                                                           %nuevo vector de variables no b�sicas // vector columna
x = x0(indb);
end
%fin de b�squeda de �ndices en caso degenerado



%c�lculo de los costos reducidos
cb = c(indb);
cnb = c(indnb);
B = A(:,indb);

PI = linsolve(B',cb);

tamnb = length(indnb);
vecost = zeros(tamnb, 1);

for i = 1: tamnb
    vecost(i) = cnb(i) - (PI')*A(:,indnb(i));                          %vector columna
end
%fin del c�lculo de costos reducidos

if vecost>= 0
    sol = 1;
    break                                                                   %%Soluci�n �ptima encontrada, rompe el while
end

% si el while no para hay algun negativo, por lo tanto el m�nimo es negativo

jestrella = indnb(1);
for i= 1:tamnb
    if vecost(i)<0
    jestrella = indnb(i);
    break
    end
end

w = linsolve(B, A(:, jestrella));
w = w(:);

if w <= 0
    sol = 2;% no hay soluci�n rompe el while
    break
end

mini = 100000000;
jp = 1;

for i = 1: length(w)
    if w(i) ~= 0
        if x(i)/w(i) < mini 
            teta = x(i)/w(i);
            jp = i;  %este es un �ndice de las b�sicas
        end
    end
end

indjp = indb(jp)
%Actualizaci�n de x0

%teta = x(jp)/w(jp)

% for i = 1: length(indb)
%     x0(indb(i)) = x0(indb(i)) - teta*w(i);
% end
% 
% x0(jestrella) = teta;                                                        %vuelves a j estrella b�sica y a jp no basica
% x0 = x0(:);                                                                  %el resto (indb, B) se actualiza al volver a entrar al while
% 
% x0

for i = 1: length(indb)
     x(i) = x(i) - teta*w(i);
end



x(jp) = teta;

indb(jp) = jestrella;

aux = find(indnb==jestrella);
indnb(aux) = jp ;
indb'
indnb'
end
x;
x0 = zeros(n,1);
vecost


if sol == 1
for i = 1: length(indb)
    x0(indb(i)) = x(i);
end
    fprintf('hubo soluci�n')
    fx =(c')*x0;
end
if sol == 2
    fprintf('Problema no acotado')
    fx = 0;
    x0 = zeros(n,1);
end
end

