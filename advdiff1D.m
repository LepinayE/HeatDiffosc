%{
Solving the 1D advection-diffusion equation for a interseasonal injection
extraction system 
T_t + USin(wt) T_x = K T_xx (non dim)

with boundary conditions:

Injection temp kept constant: 
T(0) = 1  for 0 < t < pi/w

Purely advective flux at origin during outflow:
T_x  = 0 at x=0 for pi/w < t< 2pi /w
%}

function [x,T,ti] = advdiff1D(L,U) % Extract the first solution component as T. 
    
  for i=1:10

    % Mesh over which PDE is solved
    w = 1;
    L = 10;
    ai = (i-1)*pi/w;
    bi = i *pi/w
    ti =  linspace(ai,bi,10); % time 0<t< pi/w
    res = 50;
    x = linspace(0,L,res);
    
    
    m = 0; %symmetry none as injection in LHS
    sol = pdepe(m,@pdex1pde,@pdex1ic,@bcfun,x,ti); % PDE solver
    %{
     pdex1pde = PDE equation
     pdex1ic = initial conditions
     bcfun = boundary conditions 
     %}
    
    T = sol(:,:,1);
    
    figure(1)
    surf(x,ti,T)
    hold on
    title('Temperature of fluid during cycles')
    xlabel('Distance x')
    ylabel('Time t')
    
    
    figure (2)
    hold on
    plot(x,T,'k-')
    xlabel ('x')
    ylabel('Temp')
    
    
    figure (3)
    hold on 
    plot(ti,T)
    xlabel ('time')
    ylabel('Temp')

    
    
    
    
    
  end
        % PDE
        function [c,f,s] = pdex1pde(x,t,T,dTdx)
        
        k = 10; %k needs to be found
        U = 1; % velocity u might not be constant
            
        c = 1;
        f = (k*dTdx);
        s = -U*sin(w*t)*dTdx; 
        end
    
        function T0 = pdex1ic(x) %initial condition
       
        if x>0
        T0 = 0.0001;
        else
            T0=1;
        end
        end
    
        function [pL,qL,pR,qR] = bcfun(xL,TL,xR,TR,t)
        % BCs at _L left boundary (injection) and _R right boundary (far field)
        % BCs take form p(x,t,T) + q(x,t)f(x,t,T,dTdx)= 0
       
            
            if mod(i,2) == 1
    
                pL =TL-1; 
                qL = 0;
                pR = TR-0.0001;
                qR = 0;
    
            elseif  mod(i,2) == 0 
          
                k = 10
                pL = 0; 
                qL = 1/k;
                pR = TR-0.0001;
                qR = 0;
    
            end
    
        end
   

end