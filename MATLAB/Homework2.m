clear all

% 1/4*pi*epsilon =1 

% straight line charge
% potential of straight line charge
dV = @(x,y,x1) 2.*x1./sqrt((x-x1).^2 +y.^2);
xb = linspace(-8, 8, 100);
yb = linspace(-8, 8, 100);
[x_range,y_range]=meshgrid(xb,yb);

V = [];
for i = [1:length(x_range)]
    for j = [1:length(y_range)]
        a = x_range(i,j);
        b = y_range(i,j);
        V(i,j) = integral(@(tau) dV(a, b, tau), 0, 1);
    end
end
% E field of straight line charge
[Ex,Ey]=gradient(-V);

% plots of V of straight line charge
figure(1)
contour(x_range, y_range,V)
colorbar()
xlabel('x')
ylabel('y')
axis([-2 3 -3 3])
title('Potential and E of straight line charge')
hold on
quiver(x_range,y_range,Ex,Ey,'Color','black')
hold off

figure(2)
surf(x_range,y_range,V)
colorbar()
xlabel('x')
ylabel('y')
title('Potential of straight line charge')

% plot of E of straight line charge
figure(3)
quiver(x_range,y_range,Ex,Ey,'Color','black')
axis([-2 3 -3 3])
xlabel('x')
ylabel('y')
title('E field of straight line charge')

% right angle charge
% potential of right angle charge

dVx = @(x,y,x1ra) (x1ra.^2)./sqrt((x-x1ra).^2 + y.^2 );
dVy = @(x,y,y1ra) (y1ra)./ sqrt(x.^2+(y-y1ra).^2);
V_rightangle=[];
for i=1:length(x_range)
    for j=1:length(y_range)
        a1 = x_range(i,j);
        b1 = y_range(i,j);
        V_rightangle(i,j) = integral(@(tau) dVx(a1, b1, tau), 0, 1) + integral(@(tau) dVy(a1, b1, tau), 1, 2);
    end
end
% E field of right angle charge
[Ex_ra,Ey_ra]=gradient(-V_rightangle);
        
% plots of V of rightangle charge
figure(4)
contour(x_range, y_range, V_rightangle);
title("Potential and E of right angle charge")
colorbar()
xlabel('x')
ylabel('y')
axis([-3 3 -2 4])
hold on
quiver(x_range,y_range,Ex_ra,Ey_ra,'Color','black')
hold off

figure(5)
surf(x_range,y_range,V_rightangle)
xlabel('x')
ylabel('y')
title('Potential of right angle charge')
colorbar()


% plot of E of right angle charge
figure(6)
quiver(x_range,y_range,Ex_ra,Ey_ra,'Color','black')
axis([-3 3 -2 4])
xlabel('x')
ylabel('y')
title('E field of right angle charge')


% disc charge
% potential of disc charge
dV_prime = @(x,y,z,x1,y1) x1./sqrt((x-x1).^2 +(y-y1).^2 +z.^2);
dV = @(x,y,z) integral2(@(x1,y1) dV_prime(x,y,z,x1,y1),-2,2,@(x1) -sqrt(4-x1.^2),@(x1) sqrt(4-x1.^2));
V_disc=[];
for i = 1:length(x_range)
    for j = 1:length(y_range)
        a3 = x_range(i,j);
        b3 = y_range(i,j);
        V_disc(i,j)=dV(a3,b3,0);
    end
end
% E field of disc charge
[Ex_disc,Ey_disc]=gradient(-V_disc);

% Plot of potential of disc charge
figure(7)
contour(x_range,y_range,V_disc)
colorbar()
xlabel('x')
ylabel('y')
title("Potential and E field of disc charge")
axis([-6 6 -6 6])
hold on
quiver(x_range,y_range,Ex_disc,Ey_disc)
hold off

figure(8)
surf(x_range,y_range,V_disc)
title("Potential of disc charge")
colorbar()
xlabel('x')
ylabel('y')

% plot of E of disc charge
figure(9)
quiver(x_range,y_range,Ex_disc,Ey_disc,'Color','black')
axis([-6 6 -6 6])
xlabel('x')
ylabel('y')
title('E field of disc charge')



        
