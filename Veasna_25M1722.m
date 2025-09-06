%% Define our Functions
f1 = @(x) exp(x)-x-2;
df1 = @(x) exp(x)-1;

f2 = @(x) x^5 - 8*x^4 + 44*x^3 - 91*x^2 + 85*x - 26;
df2 = @(x) 5*x^4 - 32*x^3 + 142*x^2 - 182*x + 85;

tol = 1e-6; % tolerance

%% Run for f1 and f2
bisectionMethod(f1, -4, 1, tol, "D:\M.TECH\CMTFE\Assignment2\Bisection_f1.txt");
bisectionMethod(f2, -4, 1, tol, "D:\M.TECH\CMTFE\Assignment2\Bisection_f2.txt");
newtonMethod(f1,df1,1,tol,"D:\M.TECH\CMTFE\Assignment2\Newton_f1.txt");
newtonMethod(f2,df2,1,tol,"D:\M.TECH\CMTFE\Assignment2\Newton_f2.txt");

secantMethod(f1,-4,1,tol,"D:\M.TECH\CMTFE\Assignment2\Secant_f1.txt");
secantMethod(f2,-4,1,tol,"D:\M.TECH\CMTFE\Assignment2\Secant_f2.txt");

%% ---Bisection Method---
function root = bisectionMethod(f,a0,b0,tol,data_file_name)
    if f(a0)*f(b0) >= 0
        fprintf('No root found in given interval\n');
        return;
    end  
    
    a=a0;
    b=b0;
    iter = 0;
    data = []; % store iteration and error
    while true
        
        c = (a+b)/2; % mid point
        err =abs(b-c);
        iter = iter + 1;
        % Save iteration data
        data(iter,:) = [iter, err];
        if err <=tol
        root=c;
        fprintf('The Root by Bisection is x=%.5f after %d iteration\n',c,iter);
        break;
        else
            if f(a)*f(c)<=0
                b=c;
            else
                a=c;
            end
        end
        
    end  
    
    % Write iteration data into file
    T = array2table(data,'VariableNames',{'Iter', 'Error'});
  
    writetable( T, data_file_name,'Delimiter','\t');   
end



%% ----Newton Method----

function root= newtonMethod(f,df,x0,tol,data_file_name)

   
   n=0;
   x_old =x0;
   data = []; % Initialize data to store
  
   while true
       %%Applicability Condition
     if df(x_old)==0
       error('Method Fail');
       break;
     end
       x_new = x_old-(f(x_old)/df(x_old));
       err = abs(x_new-x_old);
       %Save data in array Update Iteration
       n=n+1   ; 
       data(n,:) = [n, err];
       %Convergence check
       if(err<=tol)
           fprintf('The root by Newton is x= %.5f after %d iteration\n',x_new,n);
           root=x_new;
           break; % Stop loop
      
       else   
       x_old =x_new;
       
       end
     
    

   end
   
   T=array2table(data,'VariableNames',{'n','Error'});
   writetable(T,data_file_name,'Delimiter','\t');
   
    
end



%% ---Secant Method ---

function root = secantMethod(f,x0,x1,tol,data_file_name)

   n=0;
   x_lh = x0; % Left Hand
   x_rh = x1; % Right Hand
   data = [];
   while true
       
       if f(x_rh) == f(x_lh)
            error('Secant method failed: denominator zero');
       end
       
       x_new = x_rh -f(x_rh)*(x_rh-x_lh)/(f(x_rh)-f(x_lh));
       err = abs(x_new-x_rh);
       n=n+1;
       data(n,:)=[n,err];
       %Convergence check
       if(err<=tol)
           fprintf('The root by Secant Method is x= %.5f after %d iteration \n',x_new,n);
           root =x_new;
           break;
       else
           x_lh = x_rh; % Old  X_(k) assign to New X_(k-1);
           x_rh = x_new; % New X_(k+1) assign to New X_(n)
       end
       
       
        T=array2table(data,'VariableNames',{'n','Error'});
      writetable(T,data_file_name,'Delimiter','\t');
     
   end
   
  
end



