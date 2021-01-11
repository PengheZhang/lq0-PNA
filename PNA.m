function Out = PNA (X, Y, lam, ItMax, pars)

% This code aims at solving the sparse multinomial logistic regression with form
%
%    min_{ l(w) + lamda*\|w\|_{q,0} }, 
%
% where l(w) is the average multinomial logistic loss,
%       \|w\|_{q,0} is the nonzero groups of parameter w,
%       lamda is penalty coefficient.  

t    = tic;
[n,p] = size(X);
[~,m] = size(Y);
Yt    = Y';
m_bar = m - 1;
Y_bar = Y(:,1:m_bar);
Xt    = X';


if  nargin>4 && isfield(pars,'mu')   
    alpha = pars.mu;                 
else
    alpha = 1e0;                         
end
                          
tol = 1e-5;
% tol = 1e-3;
eps   = 0.7;
% eps = 1e-5;
sigma = 0.1;
count1 = 0;
count2 = 0;

w  = zeros( m_bar * p, 1 );                 
W  = zeros( p, m_bar );
PW = ones(n,m_bar) / m ;                
G = Xt * ( PW - Y_bar ) / n;
g = reshape( G,[ m_bar*p, 1 ] );
D = -alpha * G;

class = 1:m_bar;
I     = 1:p;
I     = I';
Iw    = 1 : m_bar*p;
Iw    = Iw';
T0    = [];
flag  = 1;
flagT = 0;

XW    = zeros(n,m_bar);                     
obj   = multi_logistic_fun(XW,Y);                          
Error = Inf;
w_new = zeros(p*m_bar,1);

fprintf('Start to run PNA......\n'); 
fprintf('Iter      Error      loss_fun \n'); 
fprintf('------------------------------\n');

for iter = 0:ItMax
    if flag
        a = sum((W+D).*(W+D),2);
        [T,~,~] = find( a>2*lam*alpha );      
    else                                                      
        T = T0;                                            
    end
    
    S = T + p*(class-1);        
    S = reshape(S,[],1);   
    Sc = setdiff(Iw,S);                
    Tc = setdiff(I,T);

    F_n = norm([g(S);w(Sc)]);
    gam = eps * F_n;
        
    if iter                                
        Error  = F_n;               
        fprintf('%3d     %5.2e     %5.2e\n',iter,Error,obj);
        if Error<tol; break; end
    end

    U = partial_hessian( PW, m_bar, n );
    
    XT1  =  Xt(T,:);
    XT2  =  X(:,T);
    nT   =  nnz(T);
    
    Ainput.U   = U;
    Ainput.m   = m_bar;    
    Ainput.gam = gam;
    Ainput.XT1 = XT1;
    Ainput.XT2 = XT2;
    Ainput.nT = nT;
    
    if flag==0 || isempty(setdiff(T,T0))    
        direct = - pcgn( 'Matvecn', Ainput, g(S), [], 1e-8);
    	w_new(S) = w(S) + direct;  
        w_new(Sc) = 0;
        
        flagT = 1;                    
    else
        if isempty(T0)
            S0 = [];
        else
            S0 = T0 + p*(class-1);                            
            S0 = reshape(S0,[],1);
        end
        
        TTc = intersect(T0,Tc);
        SSc = intersect(S0,Sc);

        h1 = hessian_w( X, Xt, U, w(SSc), T, TTc, m_bar );
        direct = - pcgn('Matvecn', Ainput, g(S) - h1, [], 1e-8);
        w_new(S) = w(S) + direct; 
        w_new(Sc) = 0;
        
    end
    
    W_new = reshape( w_new, [p, m_bar] );   
    XW_new = X(:,T)*W_new(T,:);
    
    fz_new = multi_logistic_fun(XW_new,Y);  
    
     if  obj  < fz_new - sigma  && flagT==1      
         eps = 4*eps;
         continue;
     else  
         if obj < fz_new
            count1 = count1 + 1;
         elseif obj < fz_new + 1e-5
             count2 = count2 + 1;
         end
         obj   = fz_new; 
     end
     
     if count1 >= 2
         sigma = 0.1*sigma;
         count1 = 0;
     end
     
     if count2 >=2
        break; 
     end
    
   [g,PW]    = grad_Pw(Xt,Y_bar,XW,n);
  	w        = w_new;
    W        = W_new;
    XW       = XW_new;
	T0       = T;
	
    D    = -alpha*reshape(g,[p,m_bar]);
    if isempty(T)
        WW = inf;
        DD   = D(Tc,:).*D(Tc,:);
    elseif isequal(T,I)
        WW   = W(T,:).*W(T,:);
        DD   = 0;
    else
       WW   = W(T,:).*W(T,:);
    DD   = D(Tc,:).*D(Tc,:); 
    end
   
    if  min(sum(WW,2))<2*lam*alpha || max(sum(DD,2))>=2*lam*alpha          
        flag = 1;                                                             
        if (iter && mod(iter,2)==0 && Error>1/iter) 
            alpha = 0.5*alpha; 
        end         
    else
        flag = 0;
    end
   
end
 
[g,PW]        = grad_Pw(Xt,Y_bar,XW,n);
Pwm           = 1-sum(PW,2);
Pwm           = [PW,Pwm];
[~,predict]   = max(Pwm,[],2);
[real,~,~]    = find(Yt~=0);
W = reshape(w,p,m_bar);

time        = toc(t);
Out.sparsity = nnz(sum(abs(W),2)>=1e-8);
Out.error   = Error;
Out.normd   = norm(g);
Out.time    = time;
Out.iter    = iter;
Out.alpha   = alpha;
Out.cer     = sum(real~=predict)/n;
Out.f       = obj; 
Out.g       = g;
Out.w       = w;
Out.T       = T;
Out.Pw      = PW;

fprintf('------------------------------\n');
Out

end