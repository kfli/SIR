%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Find roots      %%
%% with SIR and fsolve %%
%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% 1D function - f(x)=x-2cos(x).
% Find the root x*, 
% SIR phi(x)=x-f(x)=2cos(x).
disp('1D case')
fun1D = @(x)(x-2*cos(x));
x0 = 2;

tic;y1D = fsolve(fun1D,x0);toc
disp('fsolve solution')
disp(y1D)

phi1D = @(x)(2*cos(x));
tic;U1D = SIR(phi1D,x0,0);toc
disp('SIR solution')
disp(U1D)

%% 2D function - f1(x1,x2)=x1-cos(x2),
%                f2(x1,x2)=x2-3cos(x1).
% Find roots (x1*,x2*), 
% SIR  phi1(x1,x2)=x1-f1(x1,x2)=cos(x2),
%      phi2(x1,x2)=x2-f2(x1,x2)=3cos(x1).

disp('2D case')
fun2D = @(x)[x(1)-cos(x(2)),...
    x(2)-3*cos(x(1))];
x0 = [-2,-2];

tic;y2D = fsolve(fun2D,x0);toc
disp('fsolve solution')
disp(y2D)

phi2D = @(x)[cos(x(2)),3*cos(x(1))];

tic;U2D = SIR(phi2D,x0,1);toc
disp('SIR solution')
disp(U2D)

%%     Surfplot    %%
% Z = (f1^2+f2^2)/2 %
hold on
[X,Y] = meshgrid(-5:0.65:5);
f1 = X-cos(Y); f2 = Y-3*cos(X);
Z = (f1^2+f2^2)/2;
s = surf(X,Y,Z,'FaceAlpha',0.45);
view([-37 60])
% Point-plot of initial guess x0
plot3(x0(1),x0(2),0,'Marker','o',...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
% Point-plot of solution (x1*,x2*)
plot3(U2D(1),U2D(2),((U2D(1)-cos(U2D(2)))^2+...
    (U2D(2)-3*cos(U2D(1)))^2)/2,...
    'Marker','o',...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')
hold off
print -depsc surfFig
