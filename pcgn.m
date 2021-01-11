%%*************************************************************************
%% pcg: preconditioned conjugate gradient with left (symmetric) 
%% preconditioner for solving A*x = b. 
%%
%% b = rhs vector.
%% resnrm = norm of residual vector b-Ax. 
%%*************************************************************************

   function  [x,resnrm,solve_ok] = pcgn(Matvecname,Ainput,b,L,tol,maxit,x0) 

   [M,N] = size(b); 
   Fnorm = @(x) norm(x,'fro');
   if ~exist('maxit'); maxit = max(50,sqrt(length(b))); end
   if ~exist('tol'); tol = 1e-6; end 
   if ~exist('x0'); x0 = sparse(M,N); end
   if ~exist('L'); L = []; end
   if isempty(L); L.precond = 0; end 
   
   reltol = tol*max(1,norm(b)); 
   solve_ok = 1;    
   x = x0; 
   if (Fnorm(x) > 1e-15)
      Ap = feval(Matvecname,Ainput,x);   
      r = b-Ap;  
   else
      r = b; 
   end
   err = Fnorm(r); resnrm(1) = err; minres = err; 
%%
   z = precondfun(L,r); 
   p = z; 
   rz_old = sum(sum(r.*z)); 
%%      
%% main loop
%%
    for iter = 1:maxit 
       Ap = feval(Matvecname,Ainput,p); 
       pAp = sum(p.*Ap); 
       alpha = rz_old/pAp; 
       r = r - alpha*Ap;
       err = Fnorm(r); resnrm(iter+1) = err; 
       if (err < minres); minres = err; end
       if (err < reltol && iter >=3); break; end  
       if (rem(iter,100) == 0); fprintf('.'); end
       if (iter > 1000) 
          %%terminate when stagnation happen
          ratio = resnrm(iter-9:iter+1)./resnrm(iter-10:iter); 
          if (min(ratio) > 0.9995) && (max(ratio) < 1.0005)
             solve_ok = -1; 
             break;
          end
       end       
       x = x + alpha*p; 
       z = precondfun(L,r); 
       rz = sum(r.*z); 
       beta = rz/rz_old; 
       p = z + beta*p; 
       rz_old = rz; 
    end
%%*************************************************************************
%%*************************************************************************
   function  q = precondfun(L,r)

   precond = L.precond; 
   if (precond == 0)
      q = r; 
   else      
      error('not implemented yet'); 
   end
%%*************************************************************************