function lambda = prop2lam( X, Y, proportion, W0 )
%Given the proportion of significant features, estimate the corresponding lambda. 
[ n, p ] = size(X);
[ ~, m ] = size(Y);
prop_def = 0.1;

if p*proportion < 1
    fprintf('The proportion is too small. Default proportion = 0.1 is adopted. \n');
   proportion = prop_def;    
end

if nargin == 3
   PW = ones(n,m-1) / m;                   
   D = -X'*( PW - Y(:,1:m-1) ) / n;
   t = sum(D.*D,2)/2;
   lambda = sort( t, 'descend' );                 
   lambda = lambda( ceil( p*proportion ) ); 
else
    XW = X * W0;
    eXW = exp(XW);
    allsum = sum( eXW,2 );
    PW = eXW( :, 1:m-1 )./ allsum;
    
    D = - X'*( PW - Y( :,1:m-1 ) ) / n;
    t = sum(D.*D,2)/2;
    lambda = sort( t, 'descend' );                 
    lambda = lambda( ceil( p*proportion ) );
end

end