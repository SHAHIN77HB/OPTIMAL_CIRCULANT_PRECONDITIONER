clear
clc
%% INPUTS  this block generates the Grcar matrix of dimension n and depth k=0.1*n
n=input('n=');
b=ones(n,1);
y=fliplr(eye(n));


%% GENERATING MATRIX FROM A FUNCTION   this block generates matrix from fourier coefficients of function f
syms t;
f=exp(i*t)+2*exp(-i*t);

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
%%

u=(1/sqrt(n))*dftmtx(n);
L=u'*c*u;
L=abs(L);
c1=u'*L*u;
%% CONDITION NUMBERS

ca=cond(a);
cca=cond(inv(c)*a);



%% DISPLYES
s1=['Condition number of " An " is=', num2str(ca)];
disp(s1)

s2=['Condition number of " |Cn|^-1 An " is=', num2str(cca)];
disp(s2)




disp('---------------------------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------------------------')





%% PLOTS

e1=real(eig(a));
plot(e1,1,'*')
title('Eigenvalues of An')
axis([-(min(e1)+1) max(e1)+1 0.9 1.1 ])
figure
e2=real(eig(inv(c)*a));
plot(e2,1,'*')
title('Eigenvalues of |Cn|^-1 An')
axis([-(min(e2)+1) max(e2)+1 0.9 1.1 ])





