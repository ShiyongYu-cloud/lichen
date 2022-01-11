function P = growthpde(para,Xmin,Xmax,Tmin,Tmax,dx,dt)
%% function for solving the parabolic PDE of lichen growth
%INPUT:
% para: vector containing model parameters [a,b,sigma]
% Xmin: lower bound of the size domain
% Xmax: upper bound of the size domain
% Tmin: lower bound of the time domain
% Tmax: upper bound of the time domain
% dx: size step
% dt: time step
%OUTPUT:
%P: discrete growth curve at the temporal and spatial grid points 
%% setting up the computing grids
X  = (Xmin:dx:Xmax)';
T  = Tmin:dt:Tmax;
nx = (Xmax-Xmin)/dx; % number of spatial grids
nt = (Tmax-Tmin)/dt; % number of temporal grids
%% build up the mass matrix
a = para(1);
b = para(2);
BETA = a./sqrt(T) + b;
sigma = para(3);
%% time stepping (backward 2nd order RK)
% specifying the initial condition (i.e. the Dirac distribution)
p0 = zeros(nx+1,1);
p0(1) = 1;
P = zeros(nx+1,nt+1);
P(:,1) = p0;
p = p0;
for j = 2:nt+1
    A = makea(nx,dx,BETA(j),sigma);
    %emax = max(real(eig(A)));
    %A = A-emax*eye(size(A));
    %G = inv(eye(size(A)) - dt*A + (dt*A)^2/2 );
    G = eye(size(A)) - dt*A;
    p = G\p;
    %p = G*p;
    P(:,j) = p;
end
P(P<0.0001) = 0; % remove trival values
P = P./repmat(sum(P),size(P,1),1); % normalize to 1
return;
%%
function A = makea(nx,dx,beta,sigma)
%% function for building up the mass matrix
%INPUT
% nx: number of size grids
% dx: size step
% beta: the drift item in the SDE
% sigma: the volatile item in the SDE
%OUTPUT
% A: a (nx+1)-by-(nx+1) square matrix in the discretized PDE dp/dt = Ap
%%
A = zeros(nx+1,nx+1);
dx2 = dx*dx;
sigma2 = sigma*sigma;
A(1,1) = -beta/(2*dx)-sigma2/(2*dx2);
A(1,2) = -beta/(2*dx)+sigma2/(2*dx2);
A(nx+1,nx) = beta/(2*dx)+sigma2/(2*dx2); 
A(nx+1,nx+1) = beta/(2*dx)-sigma2/(2*dx2);
for i = 2:nx
    A(i,i-1) = beta/(2*dx)+sigma2/(2*dx2);
    A(i,i) = -sigma2/dx2;
    A(i,i+1) = -beta/(2*dx)+sigma2/(2*dx2);
end
return;