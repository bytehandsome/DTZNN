function createfigure_trace(XMatrix1, YMatrix1, ZMatrix1)
%CREATEFIGURE(XMatrix1, YMatrix1, ZMatrix1)
%  XMATRIX1:  x 数据的矩阵
%  YMATRIX1:  y 数据的矩阵
%  ZMATRIX1:  z 数据的矩阵

%  由 MATLAB 于 11-May-2022 23:37:57 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot3 的矩阵输入创建多行
plot31 = plot3(XMatrix1,YMatrix1,ZMatrix1);
set(plot31(1),'DisplayName','预期轨迹','LineWidth',2,'LineStyle','--',...
    'Color',[0.4 0.4 0.4]);
set(plot31(2),'DisplayName','实际轨迹','LineWidth',1,...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137]);

% 创建 zlabel
zlabel('z(t)');

% 创建 ylabel
ylabel('y(t)');

% 创建 xlabel
xlabel('x(t)');

view(axes1,[-37.5 30]);
hold(axes1,'off');
% 创建 legend
legend(axes1,'show');

