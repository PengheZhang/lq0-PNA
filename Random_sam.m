function [X,Y,W_star] = Random_sam(n,m,p,type,s,ro)

switch type
    case 'Independent'  
        mu = 0.5; 
        iterval_mu = 25;

        X = rand(n,p);
        y = [];
        cnum = floor(n/m);            %sample number in each class
        for i =1:m-1 
            X( (i-1)*cnum + 1 : i*cnum , iterval_mu*(i-1) + 1 : iterval_mu*i ) = X( (i-1)*cnum + 1 : i*cnum , iterval_mu*(i-1) + 1 : iterval_mu*i ) + mu*ones(cnum,iterval_mu);
            y = [y,i * ones(1,cnum)];
        end
        X( (m-1)*cnum + 1 : n , iterval_mu*(m-1) + 1 : iterval_mu*m ) = X( (m-1)*cnum + 1 : n , iterval_mu*(m-1) + 1 : iterval_mu*m ) + mu*ones(n-(m-1)*cnum,iterval_mu);
        y = [y, m * ones(1,n-(m-1)*cnum)];
        Y = sparse(1:n,y,1,n,m);
        I = randperm(n);
        X = X(I,:);
        Y = Y(I,:);
        W_star = 0;    

    case 'Corrolated'
%         rng(6,'twister')
        W_star=zeros(p,m);
        I=randperm(p);
        for i=1:m
        w_star = zeros(p,1);
        w_star(1:s) = randn(s,1);
        w_star = w_star(I); % an n-dimensional vector with s non-zero entries
        W_star(:,i)=w_star;
        end
        % Generate measurement/design matrix
        W_star(:,m)=0;

        v     = randn(n,p);
        X     = zeros(n,p); 
        X(:,1)= randn(n,1);
        for j=1:p-1
        X(:,j+1)=ro*X(:,j)+sqrt(1-ro^2)*v(:,j);
        end
        
        XW = X * W_star;
        eXW =exp(XW);
        Pw = eXW./sum(eXW,2);
        [~,class]=max(Pw,[],2);
        Y = sparse(1:n,class,1,n,m);
        
end  

end