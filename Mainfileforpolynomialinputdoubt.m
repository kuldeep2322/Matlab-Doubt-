prompt= 'IS YOUR EQUATION A POLYNOMIAL? Y/N ' ;
z=input(prompt,'s');
if(z=='N' || z=='n')
    prompt='Enter 1 for Bisetion Method \nEnter 2 for False Position Method \nEnter 3 for Fixed Point Method \nEnter 4 for Newton Raphson Method \nEnter 5 for Secant Method \n ';
    m=input(prompt);
switch m
    case 1
        f=input('Enter f(x) ');
        a=input('Enter first starting point ');
        b=input('Enter second starting point ');
        bisection(f,a,b)
    case 2
        f=input('Enter f(x) ');
        a=input('Enter first starting point ');
        b=input('Enter second starting point ');
        falseposition(f,a,b)
    case 3
        g=input('Enter g(x) ');
        a=input('Enter one starting point ');
        fixedpoint(g,a)
    case 4
        f=input('Enter f(x) ');
        f1=input('Enter f''(x) ');
        a=input('Enter one starting point ');
        newtonraphson(f,f1,a)
    case 5
        f=input('Enter f(x) ');
        a=input('Enter first starting point ');
        b=input('Enter second starting point ');
        secant(f,a,b)
    otherwise
        disp('Choose any of the options provided to you')
end
else
    n=input('Order of polynomial ');
    prompt='Enter 1 for Muller Method \nEnter 2 for Bairstow Method \n';
    m=input(prompt);
    switch m
        case 1
            muller(n)
        case 2
            bairstow(n)
        otherwise
            disp('Choose one of the options provided ')
    end
end 

function newtonraphson(f,f1,a)
merror=input('Maximum relative approximate error you want for the method: ');
miter=input('Maximum iterations you want for the method: ');
mvalue=input('Minimum closeness of f(x) to zero: ');
ix1=a;
xm=a-(f(a)/f1(a));
ym=f(xm);
iter=1;
error(1)=100;
while (abs(ym) > abs(mvalue)) && iter < miter && error(iter)>merror
    iter=iter+1;
    a=xm;
    xm=a-(f(a)/f1(a));
    ym=f(xm);
    error(iter)=abs((xm-a)/a)*100;
end
disp(xm);

if(iter>=miter)
      disp('Function terminated due to flag 3 :');
      disp(iter);
  elseif(abs(ym) <= abs(mvalue))
      disp('Function terminated due to flag 2:');
      disp(ym);
  else
      disp('Function terminated due to flag 1:');
      disp(error(iter));
end

iterations = 1:iter;
figure;
plot(iterations, error)
  xlabel('Number of iterations');
  ylabel('Relative approximate error');

  
xrange = (ix1-(xm-ix1)):0.1:(xm+(xm-ix1));
figure;
for val=1:1:size(xrange,2)
    y(val) = f(xrange(val));
end
plot(xrange,y)
    xlabel('x');
    ylabel('f(x)');
end

function fixedpoint(g,a)
merror=input('Maximum relative approximate error you want for the method: ');
miter=input('Maximum iterations you want for the method: ');
ix1=a;
xm=g(a);

iter=1;
error(1)=100;
while iter < miter & error(iter)>merror
    iter=iter+1;
    x1=xm;
    xm=g(a);
    error(iter)=abs((xm-a)/a)*100;
end
disp(xm);
if(iter>=miter)
      disp('Function terminated due to flag 3 :');
      disp(iter);
  else
      disp('Function terminated due to flag 1:');
      disp(error(iter));
end
iterations = 1:1:iter;
figure;
plot(iterations, error)
  xlabel('Number of iterations');
  ylabel('Relative approximate error');
  
  
xrange = (ix1-(xm-ix1)):0.1:(xm+(xm-ix1));
figure;
for val=1:1:size(xrange,2)
    y(val) = g(xrange(val));
end
plot(xrange,y,xrange,xrange)
    xlabel('x');
    ylabel('g(x)');
end    


function secant(f,a,b)
merror=input('Maximum relative approximate error you want for the method: ');
miter=input('Maximum iterations you want for the method: ');
mvalue=input('Minimum closeness of f(x) to zero: ');
ix1=a;
f1=(f(b)-f(a))/(b-a);
xm=b-(f(b)/f1);
ym=f(xm);
iter=1;
error(1)=abs((b-a)/a)*100;
while (abs(ym) > abs(mvalue)) & iter < miter & error(iter)>merror
    iter=iter+1;
    b=xm;
    f1=(f(b)-f(a))/(b-a);
    xm=b-(f(b)/f1);
    ym=f(xm);
    error(iter)=abs((xm-b)/b)*100;
end
disp(xm);
if(iter>=miter)
      disp('Function terminated due to flag 3 :');
      disp(iter);
  elseif(abs(ym) <= abs(mvalue))
      disp('Function terminated due to flag 2:');
      disp(ym);
  else
      disp('Function terminated due to flag 1:');
      disp(error(iter));
end
iterations = 1:iter;
figure;
plot(iterations, error)
  xlabel('Number of iterations');
  ylabel('Relative approximate error');

  
xrange = (ix1-(xm-ix1)):0.1:(xm+(xm-ix1));
figure;
for val=1:1:size(xrange,2)
    y(val) = f(xrange(val));
end
plot(xrange,y)
    xlabel('x');
    ylabel('f(x)');
end

function muller(degree)
arr=input('Enter the coeff in Decreasing Order in []:');
k0=input('Enter First Guess :');
k1=input('Enter Second Guess :');
k2=input('Enter Third Guess :');
maxits=input('Enter Maximum no.of Iterations :');
errmax=input('Enter the Maximum Relative approximate error in(%) :');
tolmax=input('Enter the tolerance value(how much f(x) is close to zero)');
strEqn=poly2sym(arr);
strVar = 'x';
F = vectorize(inline(strEqn,strVar));
d1 = k1 - k0;
d2 = k2 - k1;
del1 = (F(k1)-F(k0))./d1;
del2 = (F(k2)-F(k1))./d2;
d = (del2-del1)./(d2+d1);
i = 3;
val = [];
err(1)=100;
err(2)=100;
err(3)=100;
tol=100;

 while ((err(i) > errmax) && (i<maxits) && (tol>tolmax))
    
    b = del2+d2.*d;
    D = sqrt(b.^2-4.*F(k2).*d); 
    
    if abs(b - D) < abs(b + D)
        E = b + D;
    else
        E = b - D;
    end
    
    h = -2.*F(k2)./E;
    p = k2 + h;
    
    str = sprintf('%3u:  % 6.6f + %6.6fi % 6.6f + %6.6fi',i,real(p),imag(p),real(F(p)),imag(F(p)));
    
    tol=abs(F(p));
    val = p;
    
    i = i + 1;
    err(i)=abs((p-k2)/k2);
    
    k0 = k1;
    k1 = k2;
    k2 = p;
    d1 = k1 - k0;
    d2 = k2 - k1;
    del1 = (F(k1)-F(k0))./d1;
    del2 = (F(k2)-F(k1))./d2;
    d = (del2-del1)./(d2+d1);
    
end

if isempty(val)
    disp('The procedure was unsuccessful.')    
end

   disp(val);
   if(i>=maxits)
      disp('Program stopped due to flag 3 :');
      disp(i);
  elseif(tol<=tolmax)
      disp('Program stopped due to flag 2:');
      disp(tol);
  else
      disp('Program stopped due to flag 1:');
      disp(err(i));
  end 
  x = val-100:1:val+100;
  plot(x,F(x));
  grid on;xlabel('x');ylabel('F(x)');title('f(x)vs x');
  ax=gca;
  ax.XAxisLocation = 'origin';
  ax.YAxisLocation = 'origin';
  figure;
  i=1:1:i;
  plot(i,err);
  grid on;xlabel('i');ylabel('relative approx error');title('error vs i');
  ax=gca;
  ax.XAxisLocation = 'origin';
  ax.YAxisLocation = 'origin';
end

function falseposition(f,a,b)
merror=input('Maximum relative approximate error you want for the method: ');
miter=input('Maximum iterations you want for the method: ');
mvalue=input('Minimum closeness of f(x) to zero: ');

xrange = (a-(b-a)):0.1:(b+(b-a));
figure;
for val=1:1:size(xrange,2)
    y(val) = f(xrange(val));
end
plot(xrange,y)
    xlabel('x');
    ylabel('f(x)');

xm=a-((b-a)/(f(b)-f(a)))*f(a);
ym=f(xm);
iter=1;
error(1)=100;
pxm=xm;
while (abs(ym) > abs(mvalue)) & iter < miter & error(iter)>merror
    iter=iter+1;
    
    if(f(a)*ym<0)
        b=xm;
    else
        a=xm;
    end
    pxm=xm;
    xm=a-((b-a)/(f(b)-f(a)))*f(a);
    ym=f(xm);
    error(iter)=abs((xm-pxm)/pxm)*100;
end
disp(xm);
if(iter>=miter)
      disp('Function terminated due to flag 3 :');
      disp(iter);
  elseif(abs(ym) <= abs(mvalue))
      disp('Function terminated due to flag 2:');
      disp(ym);
  else
      disp('Function terminated due to flag 1:');
      disp(error(iter));
end
iterations = 1:iter;
figure;
plot(iterations, error)
  xlabel('Number of iterations');
  ylabel('Relative approximate error');
end

function bairstow(degree)
    coff=input('Enter the coeff in Increasing Order in []:');
    r=input('Enter First Guess :');
    s=input('Enter Second Guess :');
    maxits=input('Enter Maximum no.of Iterations :');
    errmax=input('Enter the Maximum Relative approximate error in(%) :');
    tolmax=input('Enter the tolerance value(how much f(x) is close to zero)');
    
	scoff=coff;
	sdegree=degree;
    arrc=[];
    d=degree;
    arrb=[];
    flag=0;
    count=0;
    
    while(degree>0)
    
    for j=1:maxits
	
    d=degree;
    arrb(d+1)=coff(d+1);
    arrb(d)=coff(d)+r*arrb(d+1);
    d=d-2;
    
    while(d>=0)
          arrb(d+1)=coff(d+1) + r*arrb(d+2) + s*arrb(d+3);
          d=d-1; 
    end
    
    d=degree;
    arrc(d+1)=arrb(d+1);
    arrc(d)=arrb(d)+ r*arrc(d+1);    
    d=d-2;
    
     while(d>=0)
          arrc(d+1)=arrb(d+1) + r*arrc(d+2) + s*arrc(d+3);
          d=d-1; 
     end
     
     ds = ( arrb(1)*arrc(3) - arrb(2)*arrc(2) )/( arrc(4)*arrc(2) - arrc(3)*arrc(3) );
     dr = ( arrb(1)*arrc(4) - arrb(2)*arrc(3) )/( arrc(3)*arrc(3) - arrc(2)*arrc(4) );
     r=r+dr;
     s=s+ds;
     
     err_r = abs(dr/r)*100 ;
     err_s = abs(ds/s)*100 ;
     
     if ( errmax > err_s ||  errmax > err_r )
         
         root1 = ( r + sqrt(r*r + 4*s) )/2;
         root2 = ( r - sqrt(r*r + 4*s) )/2;
         fprintf('Solution is ');
         
         if(count==0)
             fprintf(' %f %f ',root1,root2);
         end
         count=1;
         flag=1;
         break; 
     end
	 
    end
	
    if flag==0
        break;
    end
	
     for i=1:degree-1
         coff(i)=arrb(i+2);
     end
	 
     degree=degree-2;
	 
     if degree == 2
         root1 = ( -(arrb(4)) + sqrt(arrb(4)*arrb(4) - 4*arrb(3)*arrb(5)) )/(2*arrb(5));
         root2 = ( -(arrb(4)) - sqrt( arrb(4)*arrb(4) - 4*arrb(3)*arrb(5) ))/(2*arrb(5));
         fprintf('%f %f\n ',root1,root2);
         break;
		 
     elseif degree == 1
         root = -arrb(3)/arrb(4);
         fprintf('%f\n ',root);
         break;
		 
     end      
         
    end 
	
	f=zeros(1,101);
   
	
    for k=-50:50
	
        for i=1:sdegree+1
	
            f(k+51)=f(k+51)+scoff(i)*((k)^(i-1));
	 
        end
	   
    end
	   
	 x=-50:50;  
	plot(x,f);grid on;xlabel('x');ylabel('f(x)');title('f(x)vs x');
    ax=gca;
   ax.XAxisLocation = 'origin';
   ax.YAxisLocation = 'origin';
end

function bisection(f,a,b)
merror=input('Maximum relative approximate error you want for the method: ');
miter=input('Maximum iterations you want for the method: ');
mvalue=input('Minimum closeness of f(x) to zero: ');

xrange = (a-(b-a)):0.1:(b+(b-a));
figure;
for val=1:1:size(xrange,2)
    y(val) = f(xrange(val));
end
plot(xrange,y)
    xlabel('x');
    ylabel('f(x)');

xm=(a+b)/2;
ym=f(xm);
iter=1;
error(1)=100;
pxm=xm;

while (abs(ym) > abs(mvalue)) & iter < miter & error(iter)>merror
    iter=iter+1;
    if(f(a)*ym<0)
        b=xm;
    else
        a=xm;
    end
    pxm=xm;
    xm=(a+b)/2;
    ym=f(xm);
    error(iter)=abs((xm-pxm)/pxm)*100;
end
disp(xm);

if(iter>=miter)
      disp('Function terminated due to flag 3 :');
      disp(iter);
  elseif(abs(ym) <= abs(mvalue))
      disp('Function terminated due to flag 2:');
      disp(ym);
  else
      disp('Function terminated due to flag 1:');
      disp(error(iter));
end

iterations = 1:1:iter;
figure;
plot(iterations, error)
  xlabel('Number of iterations');
  ylabel('Relative approximate error');
end
