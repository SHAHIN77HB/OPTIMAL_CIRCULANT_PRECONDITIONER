function  OPTIMAL_CIRCULANT_PRECONDITIONER_function (A_n)
%% CALCULATING OPTIMAL CIRCULANT PRECONDITIONER FOR GIVEN MATRIX A_n
n=length(A_n);
 for i=1:n
     if i==1
         C(i)=A_n(1,1);
     else
         C(i)=(sum(diag(A_n,i-1))+sum(diag(A_n,i-(n+1))))/n;
     end
 end
c1=C(2:n);
c1=fliplr(c1);
c1=[C(1) c1];
C=toeplitz(c1,C);
disp('OPTIMAL CIRCULANT PRECONDITIONER IS =')
disp(C)