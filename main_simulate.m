clc; close all; clear all;
test  = 1;
type  = {'Independent','Corrolated'};                       
ro = 1/sqrt(2);   
trail = 10;

%% fix m and p, n $\in$ n0  
p     = 10000;
n0    = 1500:500:4000;                                            
m = 3;
s_ub = 1000;
s = floor(s_ub * rand);

res_n = zeros( nnz(n0),4,trail );

for i = 1:nnz(n0)
    n = n0(i);     

    for j =1:trail
        [X,Y,W_star] = Random_sam(n,m,p,type{test},s,ro);
        lamda = prop2lam( X, Y, 0.007 );                     

        out = PNA(X,Y,lamda,1000);
        res_n (i,1,j) = out.cer;
        res_n (i,2,j) = out.time;
        res_n (i,3,j) = out.f;
        res_n (i,4,j) = out.sparsity;
        
    end
end
res_n = sum(res_n,3)/trail;


 %% fix m,n  p $\in$ p0
p0   = 5000:5000:30000;
n    = 1500;                                            
m = 3;
s_ub = 1000;
s = floor(s_ub * rand);

res_p = zeros( nnz(p0),4,trail );

for i = 1:nnz(p0)
    p = p0(i);     

    for j =1:trail
        [X,Y,W_star] = Random_sam(n,m,p,type{test},s,ro);   
        lamda = prop2lam( X, Y, 0.07*s_ub/p ); 
        
        out = PNA(X,Y,lamda,1000);
        res_p (i,1,j) = out.cer;
        res_p (i,2,j) = out.time;
        res_p (i,3,j) = out.f;
        res_p (i,4,j) = out.sparsity;
        
    end
end
res_p = sum(res_p,3)/trail;
 
%% fix p,n  m $\in$ m0
m0   = 4:4:24;
n    = 1500;                                            
p    = 10000;
s_ub = 1000;
s = floor(s_ub * rand);

res_m = zeros( nnz(m0),4,trail );

for i = 1:nnz(m0)
    m = m0(i);     

    for j =1:trail
        [X,Y,W_star] = Random_sam(n,m,p,type{test},s,ro);
        lamda = prop2lam( X, Y, 0.007 );  
        
        out = PNA(X,Y,lamda,1000);
        res_m (i,1,j) = out.cer;
        res_m (i,2,j) = out.time;
        res_m (i,3,j) = out.f;
        res_m (i,4,j) = out.sparsity;
        
    end
end
res_m = sum(res_m,3)/trail

%% fix n,p,m ro $\in$ ro0 

test  = 2;
type  = {'Independent','Corrolated'};                      
trail = 10;

n = 1500; p = 10000; m = 3;
ro0 = 0.2:0.1:0.7;
s_ub = 1000;
s = floor(s_ub * rand);

res_ro = zeros(numel(ro0),4);                                                          

for i=1:numel(ro0) 
    ro = ro0(i);
    for j=1:trail
        [X,Y,W_star] = Random_sam(n,m,p,type{test},s,ro);
        lamda = prop2lam( X, Y, 0.04 );
        
        out = PNA(X,Y,lamda,1000);
        res_ro(i,1) = res_ro(i,1) + out.cer;
        res_ro(i,2) = res_ro(i,2) + out.time;
        res_ro(i,3) = res_ro(i,3) + out.f;
        res_ro(i,4) = res_ro(i,4) + out.sparsity;
    end
end

res_ro = res_ro/trail;
