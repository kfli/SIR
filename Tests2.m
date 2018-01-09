%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Find roots      %%
%% with SIR and fsolve %%
%%%%%%%%%%%%%%%%%%%%%%%%%

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
tic;U1D = SIR2(phi1D,x0,0,0.95);toc
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

tic;U2D = SIR2(phi2D,x0,1);toc
disp('SIR solution')
disp(U2D)

%% Convergence diagrams
convdiag = zeros(41,41);
phi2D = @(x)[cos(x(2)),3*cos(x(1))];
k = 0;
for i = -5:0.25:5
    k = k+1;
    l = 0;
    for j = -5:0.25:5
        l = l+1;
        x0 = [i,j];
        [U2D,iter] = SIR2(phi2D,x0,1,0.9999);
        convdiag(k,l)=iter;
    end
end

%%     Surfplot    %%
[X,Y] = meshgrid(-5:0.25:5);
f1 = X-cos(Y); f2 = Y-3*cos(X);
Z = (f1^2+f2^2)/2;

figure('DefaultAxesFontSize',16)
surf(X,Y,Z,convdiag)
zlim([-380 380])
view([-40 50])
xlabel('$$\textbf{x}_1$$','Interpreter', 'Latex')
ylabel('$$\textbf{x}_2$$','Interpreter', 'Latex')
zlabel('$$\textbf{f}{\cdot}\textbf{f/2}$$','Interpreter', 'Latex')
hcb=colorbar;
title(hcb,'iterations')
colormap(jet)
