function h = hessian_w(X,Xt,U,y,T1,T2,m)   
XT1 = Xt(T1,:);
XT2 = X(:,T2);
nT1=nnz(T1);
nT2=nnz(T2);

yr =reshape(y,nT2,m);
Xyr=XT2*yr;
h=zeros(nT1*m,1);

for i = 1:m
    B = Xyr(:,i).*U{i,i};
    for j = i+1:m
        B = B+Xyr(:,j).*U{i,j};
    end
    for k = 1:i-1
        B = B+Xyr(:,k).*U{k,i};
    end
   h((i-1)*nT1+1:i*nT1)=XT1*B;
end

end