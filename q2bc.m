function q3bc()
close all;
x1=-0.5;
y1=-0.1*sin(0.5*pi);
x2=0.5;
y2=0.1*sin(0.5*pi);
%% Square with Holes
fd=@(p) ddiff(ddiff(drectangle(p,-1,1,-1,1), ...
    dcircle(p,x1,y1,0.2)),dcircle(p,x2,y2,0.2));
pfix=[-1,-1;1,1;-1,1;1,-1];
[pts,tri]=distmesh2d(fd,@huniform,0.05,[-1,-1;1,1],pfix);


% Find boundary points with Dirichlet BCs
% dirichlet is a column vector of indices of points with Dirichlet BCs
tol=10^(-4);

dirichlet1 = find(abs(sqrt((pts(:,1)-x1).^2+(pts(:,2)-y1).^2)-0.2)<tol);
dirichlet2 = find(abs(sqrt((pts(:,1)-x2).^2+(pts(:,2)-y2).^2)-0.2)<tol);
dirichlet =[dirichlet1;dirichlet2];

%% Start numerical solution using FEM
Npts = size(pts,1);
FreeNodes = setdiff(1:Npts,dirichlet); %mesh points with unknown values of u
A = sparse(Npts,Npts);
b = sparse(Npts,1);

%% Assembly
%% The Stiffness matrix
Ntri = size(tri,1);
for j = 1:Ntri % for all triangles    
  A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:)) ...
      + stima3(pts(tri(j,:),:));
  % stima3 computes M = 0.5*|T_j|*G*G';
end

%% The Right-hand side, i.e., the load vector
for j = 1:Ntri
  b(tri(j,:)) = 0;  % for the case where f = 0
end
% Dirichlet conditions 
u = sparse(size(pts,1),1);
u(dirichlet1)=0;
u(dirichlet2)=1;
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp');
b = b - A * u;

% Computation of the solution
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

% graphic representation
figure(2)
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp');
title('solution')
colorbar;
hold on
axis ij
view(2)
end

%%
function stiff = stima3(vertices)
d = size(vertices,2);
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
Gnew = G';
for i = 1:3
   M = [1+0.5*cos(pi*vertices(i,1)) 0.5*sin(pi*vertices(i,1));...
       0.5*sin(pi*vertices(i,1)) 1+0.5*cos(pi*vertices(i,1))];
   V = exp(-(cos(2*pi*vertices(i,1)) + (vertices(i,2) ...
       - 0.1*sin(pi*vertices(i,1)))^2));
   Gnew(:,i) = V*M*Gnew(:,i);
end
  stiff = det([ones(1,d+1);vertices']) * Gnew'*G' / prod(1:d);
end



