%% Circulant Preconditioner For Analytic Functions Of Toeplitz Matrices

clear
clc
close all

%% INPUTS  this block generates the Grcar matrix of dimension n and depth k=0.1*n
n=input('Insert the matrix dimension=');
%k=input('k=');   %If you want to k be an arbitrary integer uncomment this 
k=ceil(0.1*n);       %this command sets k=0.1*n
a=gallery('grcar',n,k);
b=ones(n,1);
h=@(x) sinh(x);
y=fliplr(eye(n));

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

u=(1/sqrt(n))*dftmtx(n);
L=u'*c*u;
L=abs(L);
c1=u'*L*u;


%% CONDITION NUMBERS

ca=cond(a);
cca=cond(inv(c)*a);


cha=cond(h(a));
chcha=cond(inv(h(c))*h(a));

%% SOLVING FOR A_n

%Solve with MINRES
[xminres,flag,relres,itermin]=minres(y*a,y*b,10^-7,1000);
[m n]=lu(real(c1));
[u,flagp,relresp,iterpmin] =minres(inv(m')*y*a*inv(m),inv(m')*y*b,10^-7,1000);
xpminres=inv(m)*u;

%Solve with GMRES
[xgmres,flag,relres,itergm]=gmres(a,b,[],10^-7,size(a,1));
[xpgmres,flag,relres,iterpgm]=gmres(inv(c)*a,inv(c)*b,[],10^-7,size(a,1));

%% SOLVING FOR h(A_n)

%Solve with MINRES
[xminresh,flagh,relresh,iterhmin] =minres(y*h(a),y*b,10^-7,1000);
[m n]=lu(h(real(c1)));
[uh,flagph,relresph,iterphmin] =minres(inv(m')*y*h(a)*inv(m),inv(m')*y*b,10^-7,1000);
xpminresh=inv(m)*uh;

%Solve with GMRES
[xgmresh,flag,relres,iterhgm]=gmres(h(a),b,[],10^-7,size(a,1));
[xpgmres,flag,relres,iterphgm]=gmres(inv(h(c))*h(a),inv(h(c))*b,[],10^-7,size(a,1));

%% DISPLYES
disp('---------------------------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------------------------')


s1=['Condition number of " An " is=', num2str(ca)];
disp(s1)

s2=['Condition number of " |Cn|^-1 An " is=', num2str(cca)];
disp(s2)

%for MINRES
i11=['Number of iteration with MINRES for " An x=b " is=' , num2str(itermin)];
%for GMRES
i12=['Number of iteration with GMRES for " An x=b " is=' , num2str(itergm(2))];
disp(i11)
disp(i12)


%for MINRES
i21=['Number of iteration with MINRES for " |Cn|^-1 An x=|Cn|^-1b " is=' , num2str(iterpmin)];
%for GMRES
i22=['Number of iteration with GMRES for " |Cn|^-1 An x=|Cn|^-1b " is=' , num2str(iterpgm(2))];
disp(i21)
disp(i22)


disp('---------------------------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------------------------')

s3=['Condition number of " h(An) " is=', num2str(cha)];
disp(s3)

s4=['Condition number of " |h(Cn)|^-1 h(An) " is=', num2str(chcha)];
disp(s4)


%for MINRES
i31=['Number of iteration with MINRES for " h(An) x=b " is=' , num2str(iterhmin)];
%for GMRES
i32=['Number of iteration with GMRES for " h(An) x=b " is=' , num2str(iterhgm(2))];
disp(i31)
disp(i32)


%for MINRES
i41=['Number of iteration with MINRES for " |h(Cn)|^-1 h(An) x=|h(Cn)|^-1b " is=' , num2str(iterphmin)];
%for GMRES
i42=['Number of iteration with GMRES for " |h(Cn)|^-1 h(An) x=|h(Cn)|^-1b " is=' , num2str(iterphgm(2))];
disp(i41)
disp(i42)

%% PLOTS

e1=real(eig(a));
plot(e1,1,'*')
title('Eigenvalues of An')
%axis([-(min(e1)+1) max(e1)+1 -(min(e1)+1) max(e1)+1 ])
figure
e2=real(eig(inv(c)*a));
plot(e2,1,'*')
title('Eigenvalues of |Cn|^-1 An')
%axis([-(min(e2)+1) max(e2)+1 -(min(e2)+1) max(e2)+1 ])
figure
eh1=real(eig(h(a)));
plot(eh1,1,'*')
title('Eigenvalues of h(An)')
%axis([-(min(eh1)+1) max(eh1)+1 -(min(eh1)+1) max(eh1)+1 ])
figure
eh2=real(eig(inv(h(c))*h(a)));
plot(eh2,1,'*')
title('Eigenvalues of |h(Cn)|^-1 h(An)')
%axis([-(min(eh2)+1) max(eh2)+1 -(min(eh2)+1) max(eh2)+1 ])


