% Problem: Two-dimensional mixed convection, laminar, incompressible, steady state flow through a pipe with heated at the top and insulated at bottom.



clear all
clc
 
L = 1; % length of the pipe
H = 1; % height of the pipe
theta = 0; % inclination angle of the pipe
T_top = 100; % temperature of the heated surface
omega = 5;
 
dx = 0.2; % grid spacing in x direction
dy = 0.2; % grid spacing in y direction
x = 0:dx:L;
y = 0:dy:H;
nx = length(x);
ny = length(y);
 
u = ones(nx,ny); % x-velocity
v = ones(nx,ny); % y-velocity
p = ones(nx,ny);
rho = 1; % density
nu = 1e-3; % viscosity
g = 9.81; % acceleration due to gravity
T = ones(nx,ny); % temperature
dt = 0.01;
for i = 2:nx-1
    for j = 2:ny-1
        % x-momentum equation
        u_xx = (u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2;
        u_yy = (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2;
        v_xy = (v(i,j+1)-v(i,j-1))/(2*dy);
        pressure_x = (p(i+1,j)-p(i,j))/dx;
        Fx = ((-rho)*(u_xx+u_yy))-(pressure_x+(rho*g*sin(theta)));
        u(i,j) = u(i,j) + dt*(1/nu*Fx-v_xy);
    end
end
 
 
 
p = ones(nx,ny); % pressure
% maxiter = 1000; % maximum number of iterations
% tol = 1e-6; % tolerance
% for iter = 1:maxiter
    % solve for pressure
    for i = 2:nx-1
        for j = 2:ny-1
            p_xx = (p(i+1,j)-(2*p(i,j))+p(i-1,j))/(dx^2);
            p_yy = (p(i,j+1)-(2*p(i,j))+p(i,j-1))/(dy^2);
            b = rho*(1/dt*(u(i,j)-u(i-1,j))/dx+(v(i,j+1)-v(i,j))/dy-p_xx-p_yy);
            p(i,j) = (1-omega)*p(i,j)+omega/(2*((1/(dx^2))+(1/(dy^2))))*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-((dx^2)*b));
        end
    end
% end
 
% solve for velocity
u_old = u;
v_old = v;
for i = 2:nx-1
    for j = 2:ny-1
        u_xx = (u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2;
        u_yy = (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2;
        v_xx = (v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2;
        v_yy = (v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2;
        pressure_x = (p(i+1,j)-p(i,j))/dx;
        pressure_y = (p(i,j+1)-p(i,j))/dy;
        Fx = ((-rho)*(u_xx+u_yy))-(pressure_x+(rho*g*sin(theta)));
        Fy = ((-rho*(v_xx+v_yy))-(pressure_y));
        u(i,j) = u(i,j) + dt*((1/nu)*(Fx-v(i,j)/dx)*((u(i,j)-u(i-1,j))/dx-v(i,j)/dy)*((u(i,j)-u(i,j-1))/dy));
        v(i,j) = v(i,j) + dt*((1/nu)*(Fy-u(i,j)/dx)*((v(i,j)-v(i-1,j))/dx-v(i,j)/dy)*((v(i,j)-v(i,j-1))/dy));
    end
end
 
%  check for convergence
% if max(max(abs(u-u_old))) < tol && max(max(abs(v-v_old))) < tol
%     return;
% end
 
% no-slip condition
u(1,:) = 0;
u(end,:) = 0;
v(1,:) = 0;
v(end,:) = 0;
 
% inlet and outlet boundary conditions
u(:,1) = u(:,2);
u(:,end) = u(:,end-1);
v(:,1) = 0;
% v(:,end) = 0;
 
% heating condition at top surface
T(end,:) = T_top;
 
% calculate velocity magnitude
u_mag = sqrt((u.^2)+(v.^2));
 
% plot velocity magnitude and temperature
 
figure
contourf(x,y,u_mag,7)
colorbar
xlabel('x')
ylabel('y')
title('Velocity Magnitude')
 
figure
contourf(x,y,T)
colorbar
xlabel('x')
ylabel('y')
title('Temperature')
