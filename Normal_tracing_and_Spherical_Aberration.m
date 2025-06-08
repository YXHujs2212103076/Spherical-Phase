clear ;
clc ;
set(0, 'DefaultAxesFontName', 'Microsoft YaHei');
set(0, 'DefaultTextFontName', 'Microsoft YaHei');
% R = Lens radius (μm)
% F = Focal length (μm)
% FLAE = Focal Length Absolute Error
% 
clear;
clc;
period=0.52e-6;
height=0.75e-6;
lamda=1.1e-6;

a=2;%The coefficient of the ratio of focal length to radius

Number_LENS_D=97;%Number of structural units
R = Number_LENS_D*period/2;
x = -(Number_LENS_D-1)*period/2:period:(Number_LENS_D-1)*period/2;
y = -(Number_LENS_D-1)*period/2:period:(Number_LENS_D-1)*period/2;

% circular area
target_h=zeros(Number_LENS_D,Number_LENS_D);
for m=1:Number_LENS_D
    for n=1:Number_LENS_D
        if sqrt((x(m)-0)^2+(y(n)-0)^2)<R
            target_h(m,n)=1;
        end
    end
end

% 
target_x=zeros(Number_LENS_D,Number_LENS_D);
f_x =a*R;%=============================================================================
NAx=sin(atan(R/f_x));%Spherical phase numerical aperture
dx=1.22*lamda/(NAx);% diffraction-limited,DL

for m=1:Number_LENS_D
    for n=1:Number_LENS_D
            if target_h(m,n)>0
                 target_x(m,n)=(-(2*pi)/lamda)*(  f_x-sqrt( f_x^2-((x(m)-0)^2+(y(n)-0)^2) )   );%Spherical phase distribution
            end
    end                                                          
end
Lx =target_x((Number_LENS_D+1)/2,:)/(-(2*pi)/lamda);%Spherical phase optical path distribution, displaying the lens radius axis
Tar_x=Lx;
dP_dx = gradient(Lx,x);%

for i=1:m
    point_tangent(i)=dP_dx(1,i);%Tangent slope
    point_normal(i) = -1/point_tangent(i);%Normal slope
    dleta(i)= point_normal(i)*(0 - x(i))+Lx(1,i);%Intersection point of normal oblique equation and coordinate axis
end

b=floor (4*(Number_LENS_D+1)/5) ;%A randomly set numerical value used to represent a point on the lens
for i=1:m
    if i>48 && i<(b+1)
        x1(i-48)=x(i);
        Dx(i-48)=point_normal(b).*(x(i)- x(b))+Lx(1,b);% Normal equation of point b
    end
    if i>(b-30) && i<(b+30)% Tangent display area
        x2(i-b+30)=x(i);
        Ex(i-b+30)=point_tangent(b).*(x(i)- x(b))+Lx(1,b);% Tangent equation of point b
    end
end  

dleta2=abs(dleta-f_x);% Focal Length Absolute Error，FLAE
ratio_x=dleta2/dx;%-------------------------,FLAE/DL


target_y=zeros(Number_LENS_D,Number_LENS_D);
f_y =a*R;%=========================================================================
NAy=sin(atan(R/f_y));%Hyperbolic phase numerical aperture
dy=1.22*lamda/(NAy);%Diffraction Limit，DL

for m=1:Number_LENS_D
    for n=1:Number_LENS_D
            if target_h(m,n)>0
                target_y(m,n)=(-(2*pi)/lamda)*(sqrt((x(m)-0)^2+(y(n)-0)^2+f_y^2)-f_y);%Hyperbolic phase distribution
            end
    end
end

Ly =target_y((Number_LENS_D+1)/2,:)/(-(2*pi)/lamda);%Hyperbolic phase optical path distribution, displaying the lens radius axis
Tar_y=Ly;
dP_dy = gradient(Ly,y);
for i=1:m
    point_tangent(i)=dP_dy(1,i);%Tangent slope
    point_normal(i) = -1/point_tangent(i);%Normal slope
    dleta(i)= point_normal(i)*(0 - y(i))+Ly(1,i);%Intersection point of normal oblique equation and coordinate axis
end
for i=1:m
    if i>48 && i<(b+1)
        y1(i-48)=y(i);
        Dy(i-48)=point_normal(b).*(y(i)- y(b))+Ly(1,b);% Normal equation of point b
    end
    if i>(b-30) && i<(b+30) % Tangent display area
        y2(i-b+30)=y(i);
        Ey(i-b+30)=point_tangent(b).*(y(i)- y(b))+Ly(1,b);% Tangent equation of point b
    end
end  

dleta3=abs(dleta-f_y);% Focal Length Absolute Error，FLAE
ratio_y=dleta3/dy;%-------------------------，FLAE/DL


a1=20;%Drawing subgraph title font size
a2=20;%Draw other font sizes
set(groot, 'DefaultAxesFontSize', a2);  % 
set(groot, 'DefaultAxesFontWeight', 'bold');
set(groot, 'DefaultAxesLineWidth', 1.5); % 

figure('WindowState', 'maximized');
set(gcf, 'Color', 'w');

subplot(1,2,1)
plot(y*1e6,Tar_y*1e6 , 'Color', [0 0.7 0], 'LineWidth', 2);%
hold on
plot(  y1*1e6,Dy*1e6 , 'Color', [0 0 1], 'LineWidth', 2);%
hold on
plot(  y2*1e6,Ey*1e6 , '--','Color', [0 0 0], 'LineWidth', 2 );%
x_mark = y1(1);
y_mark = Dy(1)*1e6;
text(x_mark, y_mark, sprintf('(%0.2f, %0.2f)', x_mark, y_mark), ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right', ...
    'FontSize', a2, ...
    'Color', 'r');
x_mark = y(b)*1e6;
y_mark = Tar_y(b)*1e6;
plot(x_mark, y_mark, 'p', 'MarkerSize', a2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); 
xlabel('r (um)', 'FontSize', a2, 'Color', 'k', 'FontWeight', 'bold');
ylabel('z (um)', 'FontSize', a2, 'Color', 'k', 'FontWeight', 'bold');
title('Hyperbolic', 'FontSize', a1,  'Color', 'k', 'FontWeight', 'bold');
text(-0.1, 1.1, '(a)', 'Units', 'normalized', 'FontSize', a2,  'Color', 'k', 'FontWeight', 'bold');

subplot(1,2,2)
plot( x*1e6,Tar_x*1e6  , 'Color', [0 0.7 0], 'LineWidth', 2);
hold on
plot(  x1*1e6,Dx*1e6 , 'Color', [1 0 0], 'LineWidth', 2);
hold on
plot(  x2*1e6,Ex*1e6 , '--','Color', [0 0 0], 'LineWidth', 2 );
x_mark = x1(1);
y_mark = Dx(1)*1e6;
text(x_mark, y_mark, sprintf('(%0.2f, %0.2f)', x_mark, y_mark), ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right', ...
    'FontSize', a2, ...
    'Color', 'r');
x_mark = x(b)*1e6;
y_mark = Tar_x(b)*1e6;
plot(x_mark, y_mark, 'p', 'MarkerSize', a2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); % 
xlabel('r (um)', 'FontSize', a2, 'Color', 'k', 'FontWeight', 'bold');
ylabel('z (um)', 'FontSize', a2, 'Color', 'k', 'FontWeight', 'bold');
title('Sphere', 'FontSize', a1,  'Color', 'k', 'FontWeight', 'bold');
text(-0.1, 1.1, '(b)', 'Units', 'normalized', 'FontSize', a2,  'Color', 'k');

h = suptitle('Normal-tracing diagram of equiphase surfaces.(F=50.44um,F/R=2)');
set(h, 'FontWeight', 'bold', 'FontSize', a2);


figure('WindowState', 'maximized');
set(gcf, 'Color', 'w');

subplot(1,2,1)
plot(y*1e6,ratio_y,'-k', 'LineWidth', 2);
hold on
x_mark = y(b)*1e6;
y_mark = ratio_y(b);
plot(x_mark, y_mark, 'p', 'MarkerSize', a2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); % 
xlabel('r (um)', 'FontSize', a2, 'Color', 'k', 'FontWeight', 'bold');
ylabel('Spherical Aberration', 'FontSize', a2, 'Color', 'k', 'FontWeight', 'bold');
title('Hyperbolic', 'FontSize', a1,  'Color', 'k', 'FontWeight', 'bold');
text(-0.1, 1.1, '(a)', 'Units', 'normalized', 'FontSize', a2,  'Color', 'k');

subplot(1,2,2)
plot(x*1e6,ratio_x,'-k', 'LineWidth', 2);
hold on
x_mark = x(b)*1e6;
y_mark = ratio_x(b);
plot(x_mark, y_mark, 'p', 'MarkerSize', a2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); % 
xlabel('r (um)', 'FontSize', a2, 'Color', 'k', 'FontWeight', 'bold');
ylabel('Spherical Aberration', 'FontSize', a2, 'Color', 'k', 'FontWeight', 'bold');
title('Sphere', 'FontSize', a1,  'Color', 'k', 'FontWeight', 'bold');
text(-0.1, 1.1, '(b)', 'Units', 'normalized', 'FontSize', a2,  'Color', 'k');

h = suptitle('Spherical Aberration(F=50.44um,F/R=2)');
set(h, 'FontWeight', 'bold', 'FontSize', a2);