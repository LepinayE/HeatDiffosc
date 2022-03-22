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

%{
for i = 1:2 
    
    z = advdiff(10,1)
end
%}
i = 1

function [x,T,t] = advdiff(L,U) % Extract the first solution component as T.  
  % Mesh over which PDE is solved
        w = 2;
        i = 1
        L = 10;
        res = 50;
        x = linspace(0,L,res);

        ai = (i-1)*pi/w;
        bi = i *pi/w;
        t =  linspace(ai,bi,10);

        m = 0; %symmetry none as injection in LHS
        
       
        sol = pdepe(m,@pdex1pde,@pdex1ic,@bcfun,x,t); % PDE solver
        %{
         pdex1pde = PDE equation
         pdex1ic = initial conditions
         bcfun = boundary conditions 
         %}
        
        T = sol(:,:,1);
        AqT = T(end,:)

        d1 = size(AqT)
        d2 = size (x)
        
        figure(1)
        surf(x,t,T)
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
        plot(t,T)
        xlabel ('time')
        ylabel('Temp')
    
        figure (4)
        hold on 
        plot(x,T(end,:))
        xlabel ('x')
        ylabel('Temp')   
   

    

% PDE
function [c,f,s] = pdex1pde(x,t,T,dTdx)
  

k = 10; %k needs to be found
U = 1; % velocity u might not be constant
    
c = 1;
f = (k*dTdx);
s = -U*sin(w*t)*dTdx; 
end
    
%initial condition

function T0 = pdex1ic(x)
 if i == 1
  
        if x>0
            T0 = 0.0001;
        else
            T0=1;
        end

 elseif i > 1
      T0 = T(end,:);
 end

end
 
%BCs
function [pL,qL,pR,qR] = bcfun(xL,TL,xR,TR,t)
% BCs at _L left boundary (injection) and _R right boundary (far field)
% BCs take form p(x,t,T) + q(x,t)f(x,t,T,dTdx)= 0

    if  mod(i,2) == 1;
       
        pL =TL-1; 
        qL = 0;
        pR = TR-0.0001;
        qR = 0;
        
    else mod(i,2) == 0;
        k = 10;
        pL = 0; 
        qL = 1/k;
        pR = TR-0.0001;
        qR = 0;
    end


end
end

