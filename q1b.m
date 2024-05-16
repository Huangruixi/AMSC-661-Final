function q1b()

h = 0.001;
Tmax = 2*10^6;
Nsteps = ceil(Tmax/h);
tol = 1.0e-14;
itermax = 20;
tic % start measuring CPU time
[sol,t_vector] = DIRK2_adaptive_step(h,tol,itermax,Tmax);
toc

figure;
subplot(3,1,1);
plot(t_vector,sol(:,1));
xlabel('t');
ylabel('x');
title('miu=10^6');

subplot(3,1,2);
plot(t_vector,sol(:,2));
xlabel('t');
ylabel('y');

subplot(3,1,3);
plot(sol(:,1),sol(:,2));
xlabel('x');
ylabel('y');
end

% the right-hand side
function dy = func(y)
    miu=10^6;
    dy = zeros(2,1);
    dy(1) = y(2);
    dy(2) = miu*(1-y(1)^2)*y(2)-y(1);
end

% the Jacobian matrix for the right-hand side
function J = Jac(y)
    miu=10^6;
    J = zeros(2);
    J(1,1) = 0;
    J(1,2) = 1;
    J(2,1) = -2*miu*y(1)*y(2)-1;
    J(2,2) = miu*(1-y(1)^2);
end

%% DIRK2

function knew = NewtonIterDIRK2(y,h,k,gamma)
    aux = y + h*gamma*k;
    F = k - func(aux);
    DF = eye(2) - h*gamma*Jac(aux);
    knew = k - DF\F;
end

function [sol,t_vector] = DIRK2_adaptive_step(h,tol,itermax,Tmax)
    
    t = 0;
    t_vector = t;
    y0=[2 0];
    sol=y0;
    y = y0';
    while t+h<Tmax+h
       t = t+h;
       t_vector = [t_vector;t];
       [k1,k2,y,ynew] = DIRK2step(y,h,tol,itermax,Tmax);
       Error = norm(h*(k1 - k1/2 - k2/2));
       epsilon = (1e-5)+(1e-5)*norm(ynew);
       h = h*min(5,max(0.5,0.8*(epsilon/Error)^(1/2)));
       if Error > epsilon
           while Error > epsilon
                [k1,k2,y,ynew] = DIRK2step(y,h,tol,itermax,Tmax);
                Error = norm(h*(k1 - k1/2 - k2/2));
                epsilon = (1e-5)+(1e-5)*norm(ynew);
                h = h*min(5,max(0.5,0.8*(epsilon/Error)^(1/2)));
           end
       end
       y=ynew;
       sol = [sol;y'];
       
    end    
    end

function [k1,k2,y,ynew] = DIRK2step(y,h,tol,itermax,Tmax)
    gamma = 1.0 - 1.0/sqrt(2);
    k1 = func(y);
    for j = 1 : itermax
        k1 = NewtonIterDIRK2(y,h,k1,gamma);
        if norm(k1 - func(y + h*gamma*k1)) < tol
            break
        end
    end
    k2 = k1;
    y = y + h*(1-gamma)*k1;
    for j =1 : itermax
        k2 = NewtonIterDIRK2(y,h,k2,gamma);
        aux = y + h*gamma*k2;
        if norm(k2 - func(aux)) < tol
            break
        end
    end
    ynew = aux;
end

