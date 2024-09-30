function [simpint,msimp] =simp(a,b,fun,tol)
simpint=0;
iter=0;
msimp=1;
while abs(simpint-exp(b)+1)>=tol
h=(b-a)/msimp;
x=[a:h/2:b];
dim=length(x);
y=feval(fun,x);
if size(y)==1
    y=diag(ones(dim))*y;
end
simpint=(h/6)*(y(1)+2*sum(y(3:2:2*msimp-1))+4*sum(y(2:2:2*msimp))+y(2*msimp+1));
msimp=msimp+1;
end
return 