%% FCM Project 3
%Jonathan Engle
%Due March 12
%% Define Functions and parameters
clear;clc;clear all 
tic
format long
%lower bound 
global a
a=0;
%upper bound
global b
b=3;
% Global alpha
global alpha
alpha=1/3;
%Global tol
global tol 
%tol=0.01; %0.0001
tol=0.01;

%function
global fun 
fun=@(t)(exp(t));

%% Composite trapezpodal rule
[inttrap,mtrap] =trap(a,b,fun,tol);
error=abs(inttrap-exp(b)+1);
mtrap;
traptheory=abs(-(b-a)*(mtrap^2/(12))*exp(3));
%% Composite midpoint
[int,iter,mcomp] =compmidpoint(a,b,fun,tol);
error=abs(int-exp(b)+1);
int;
iter;
Midtheory=abs(-(b-a)*(iter^2/(24))*exp(3));
%% Composite Simpsons First rule
[simpint,msimp] =simp(a,b,fun,tol);
msimp;
Simptheory=abs(-(b-a)*(1*msimp^4/(2880))*exp(3));
%% Global refinement algorithms
three=((a+b)/2);
inthree=feval(fun,three)*(b-a);
iterthree=0;
three=0;
m=1;

% while abs(inthree-exp(3)+1)>=tol
% h=(b-a)/m;
% x=a+1/2:h:b;
% summand=0;
% for i=1:m
%     summand = summand+feval(fun,a+(i-1)*h+h/6)+feval(fun,a+(i-1)*h+5*h/6);
% end
% three=inthree;
% inthree=1/3*(inthree+h*summand);
% m=m*3;
% iterthree=iterthree+1;
% 
% end
% r=log((inthree-three)/((exp(3)-1)-inthree)+1)/log(3);
% error3=abs(inthree-exp(3)+1);
% iterthree;

three=ones(11,1);
inthree=ones(11,1);
for j=1:11
h=(b-a)/m;
x=a+1/2:h:b;
% summand=zeros(m,1);
summand=0;
for i=1:m
    summand = summand+feval(fun,a+(i-1)*h+h/6)+feval(fun,a+(i-1)*h+5*h/6);
end
three(j+1)=inthree(j);
inthree(j+1)=1/3*(inthree(j)+h*summand);
m=m*3;
iterthree=iterthree+1;
r(j)=log((inthree(j)-three(j))/((exp(3)-1)-inthree(j))+1)/log(3);
end
error3=abs(inthree-exp(3)+1);
inthree;
iterthree;
% plot(error3)

hold on
plot(log(error3))
hold off

toc