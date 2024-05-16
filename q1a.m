function q1a()
h = 0.001; % time step
Tmax = 200; 
Nsteps = ceil(Tmax/h);
tol = 1.0e-14;
itermax = 20;

y0 = [2,0.0]'; % initial condition

sol = zeros(Nsteps+1,2);
t = h*(1:(Nsteps+1))';
sol(1,:) = y0;
tic  % start measuring CPU time

method_name = "DIRK2";
for j = 1 : Nsteps
    sol(j+1,:) = DIRK2step(sol(j,:)',h,tol,itermax)';
end
toc % end measuring CPU time

figure;
fsz = 20; % fontsize
subplot(3,1,1);
plot(t,sol(:,1));
xlabel('t');
ylabel('x');
title('miu=100');

subplot(3,1,2);
plot(t,sol(:,2));
xlabel('t');
ylabel('y');

subplot(3,1,3);
plot(sol(:,1),sol(:,2));
xlabel('x');
ylabel('y');

end

% the right-hand side
function dy = func(y) 
    miu=100;
    dy = zeros(2,1);
    dy(1) = y(2);
    dy(2) = miu*(1-y(1)^2)*y(2)-y(1);
end

% the Jacobian matrix for the right-hand side
function J = Jac(y)
    miu=100;
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

function ynew = DIRK2step(y,h,tol,itermax)
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

