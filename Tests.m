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

options = optimset('Display','off','TolFun', 1e-8, 'MaxIter', 150);
tic;y1D = fsolve(fun1D,x0,options);toc
disp(['fsolve solution: ', num2str(y1D),'; Residual: ',...
    num2str(sum(fun1D(y1D).^2))])

phi1D = @(x)(2*cos(x));
tic;U1D = SIR(phi1D,x0,0);toc
disp(['SIR solution: ', num2str(U1D),'; Residual: ',...
    num2str(sum(fun1D(U1D).^2))])
fprintf(1, '\n');

%% 2D function - f1(x1,x2)=x1-cos(x2),
%                f2(x1,x2)=x2-3cos(x1).
% Find roots (x1*,x2*), 
% SIR  phi1(x1,x2)=x1-f1(x1,x2)=cos(x2),
%      phi2(x1,x2)=x2-f2(x1,x2)=3cos(x1).
disp('2D case')
fun2D = @(x)[x(1)-cos(x(2)),...
    x(2)-3*cos(x(1))];
x0 = [-2,-2];

tic;y2D = fsolve(fun2D,x0,options);toc
disp(['fsolve solution: ','[',num2str(y2D(1)),', ', num2str(y2D(2)),']','; Residual: ',...
    num2str(sum(fun2D(y2D).^2))])

phi2D = @(x)[cos(x(2)),3*cos(x(1))];
tic;U2D = SIR(phi2D,x0,1);toc
disp(['SIR solution: ','[',num2str(U2D(1)),', ', num2str(U2D(2)),']','; Residual: ',...
    num2str(sum(fun2D(U2D).^2))])


%% Convergence diagrams
convdiag = zeros(51,51);
k = 0;
for i = -5:0.2:5
    k = k+1;
    l = 0;
    for j = -5:0.2:5
        l = l+1;
        x0 = [i,j];
        [U2D,iter]  = SIR(phi2D,x0,1);
        if iter > 100
            iter = 150;
        end
        convdiag(k,l)=iter;
    end
end

%%     Surfplot    %%
[X,Y] = meshgrid(-5:0.2:5);
f1 = X-cos(Y); f2 = Y-3*cos(X);
Z = (f1.*f1+f2.*f2)/2;

figure(1)
surf(X,Y,Z,convdiag)
view([-40 50])
title('Convergence diagram: SIR-s [f_1 = x_1-cos(x_2), f_2 = x_2-3*cos(x_1)]')
xlabel('$$\textbf{x}_1$$','Interpreter', 'Latex')
ylabel('$$\textbf{x}_2$$','Interpreter', 'Latex')
zlabel('$$\textbf{f}{\cdot}\textbf{f}/2$$','Interpreter', 'Latex')
caxis([0 150])
hcb = colorbar;
hcb.Ticks = [0,50,100,150];
hcb.TickLabels = {'0 it.','50 it.','100 it.','No conv'};  
colormap(jet)

