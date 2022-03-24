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

IC = [];

for j = 1:10
    global w;
    % Mesh over which PDE is solved
    w = 1;
    L = 10;
    aj = (j-i)*pi/w;
    bj = j * pi/w ;
    t =  linspace(aj,bj,10); % time 0<t< pi/w
    res = 50;
    x = linspace(0,L,res);
    
    
    m = 0; %symmetry none as injection in LHS

    % PDE solver
    % Need to specify to pdepe that the Ic and BC functions are not
    % changing all their variables

    sol = pdepe(m,@pdex1pde, @(x) pdex1ic(x, IC), @(xL,TL,xR,TR,t) bcfun(xL,TL,xR,TR,t,j), x, t); 
    %{
     pdex1pde = PDE equation
     pdex1ic = initial conditions
     bcfun = boundary conditions 
     %}
    
    T = sol(:,:,1);
    
    newIC = T(end,:);
    
 

    subplot(1,3,1)
    surf(x,t,T)
    hold on
    title('Temperature of fluid during cycles')
    xlabel('Distance x')
    ylabel('Time t')
    
    
    subplot(1,3,2)
    hold on
    plot(x,T,'k-')
    xlabel ('x')
    ylabel('Temp')
    
    
    subplot(1,3,3)
    hold on 
    plot(t,T)
    xlabel ('time')
    ylabel('Temp')
 
    
    drawnow

end

% PDE
function [c,f,s] = pdex1pde(x,t,T,dTdx)
    global w;
k = 10; %k needs to be found
U = 1; % velocity u might not be constant
    
c = 1;
f = (k*dTdx);
s = -U*sin(w*t)*dTdx; 
end

%initial condition
function T0 = pdex1ic(x, newIC) 
    if isempty(newIC)
        if x>0
        T0 = 0.0001;
        else
            T0=1;
        end

    else
        X = linspace(0,10,25)
        T0 = interpl(X, newIC, x);
    end
end

%BCs
function [pL,qL,pR,qR] = bcfun(xL,TL,xR,TR,t,j)
% BCs at _L left boundary (injection) and _R right boundary (far field)
% BCs take form p(x,t,T) + q(x,t)f(x,t,T,dTdx)= 0

    if mod(j,2) == 1
    
        pL =TL-1; 
        qL = 0;
        pR = TR-0.0001;
        qR = 0;
    
    elseif  mod(j,2) == 0 
  
        k = 10
        pL = 0; 
        qL = 1/k;
        pR = TR-0.0001;
        qR = 0;

    end


end