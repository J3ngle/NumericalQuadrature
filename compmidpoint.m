function [int,iter,mcomp] =compmidpoint(a,b,fun,tol)
    int=0;
    iter=1;
    mcomp=2;
    while abs(int-exp(b)+1)>=tol
        h=(b-a)/iter;
        x=a+h/2:h:b;
        dim=length(x);
        y=feval(fun,x);
        if size(y)==1
            y=diag(ones(dim))*y;
        end
        int=h*sum(y);
        iter=iter+1;
    end
end