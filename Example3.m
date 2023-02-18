%% Circulant Preconditioner For Analytic Functions Of Toeplitz Matrices

clear
clc
close all
%% Inputs 

n=input('Insert the matrix dimension=');
b=ones(n,1);
h=@(x) exp(x);


%% GENERATING MATRIX FROM A FUNCTION   this block generates matrix from fourier coefficients of function f
syms t;
f=t^2;

for j=1:n
   aa(j)=(1/(2*pi))*int(f*exp(-i*t*(j-1)),-pi,pi); 
end
for k=1:n
    bb(k)=(1/(2*pi))*int(f*exp(-i*t*-(k-1)),-pi,pi);
end
aa=double(aa);
bb=double(bb);

a=toeplitz(aa,bb);

%% CALCULATING OPTIMAL PRECONDITIONER
 
 for i=1:n
     if i==1
         c(i)=a(1,1);
     else
         c(i)=(sum(diag(a,i-1))+sum(diag(a,i-(n+1))))/n;
     end
end

c1=c(2:n);
c1=fliplr(c1);
c1=[c(1) c1];

c=toeplitz(c1,c);



%% CONDITION NUMBERS

ca=cond(a);
cca=cond(inv(c)*a);


cha=cond(h(a));
chcha=cond(inv(h(c))*h(a));

%% SOLVING FOR A_n


[xpcg,flag,relres,iter] =pcg(a,b,10^-7,1000);
[m, n]=chol(c);
[u,flagp,relresp,iterp] =pcg(inv(m')*a*inv(m),inv(m')*b,10^-7,1000);
xppcg=inv(m)*u;

%% SOLVING FOR h(A_n)

[xpcgh,flagh,relresh,iterh] =pcg(h(a),b,10^-7,1000);
[m, n]=chol(h(c));
[uh,flagph,relresph,iterph] =pcg(inv(m')*h(a)*inv(m),inv(m')*b,10^-7,1000);
xppcgh=inv(m)*uh;

%% DISPLYES
disp('---------------------------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------------------------')



s1=['Condition number of " An " is=', num2str(ca)];
disp(s1)

s2=['Condition number of " |Cn|^-1 An " is=', num2str(cca)];
disp(s2)

i1=['Number of iteration with CG for " An x=b " is=' , num2str(iter)];
disp(i1)

i2=['Number of iteration with CG for " |Cn|^-1 An x=|Cn|^-1b " is=' , num2str(iterp)];
disp(i2)

disp('---------------------------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------------------------')

s3=['Condition number of " h(An) " is=', num2str(cha)];
disp(s3)

s4=['Condition number of " |h(Cn)|^-1 h(An) " is=', num2str(chcha)];
disp(s4)

i3=['Number of iteration with CG for " h(An) x=b " is=' , num2str(iterh)];
disp(i3)

i4=['Number of iteration with CG for " |h(Cn)|^-1 h(An) x=|h(Cn)|^-1b " is=' , num2str(iterph)];
disp(i4)

%% PLOTS

e1=real(eig(a));
plot(e1,1,'*')
title('Eigenvalues of An')
axis([-(min(e1)+1) max(e1)+1 -(min(e1)+1) max(e1)+1 ])
figure
e2=real(eig(inv(c)*a));
plot(e2,1,'*')
title('Eigenvalues of |Cn|^-1 An')
axis([-(min(e2)+1) max(e2)+1 -(min(e2)+1) max(e2)+1 ])
figure
eh1=real(eig(h(a)));
plot(eh1,1,'*')
title('Eigenvalues of h(An)')
axis([-(min(eh1)+1) max(eh1)+1 -(min(eh1)+1) max(eh1)+1 ])
figure
eh2=real(eig(inv(h(c))*h(a)));
plot(eh2,1,'*')
title('Eigenvalues of |h(Cn)|^-1 h(An)')
axis([-(min(eh2)+1) max(eh2)+1 -(min(eh2)+1) max(eh2)+1 ])
