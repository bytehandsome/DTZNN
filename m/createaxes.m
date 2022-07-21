function createaxes(Parent1, X1, YMatrix1)
%CREATEAXES(Parent1, X1, YMatrix1)
%  PARENT1:  axes parent
%  X1:  x 数据的向量
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 11-May-2022 15:35:23 自动生成

% 创建 axes
axes1 = axes('Parent',Parent1,...
    'Position',[0.317222222222221 0.255358565737051 0.309999999999999 0.31]);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',2);
set(plot1(1),'DisplayName','X','Color',[0 0 0]);
set(plot1(2),'DisplayName','Y','LineStyle','--','Color',[0.1 0.1 0.1]);
set(plot1(3),'DisplayName','Z','LineStyle',':','Color',[0.2 0.2 0.2]);
set(plot1(4),'DisplayName','S','LineStyle','-.','Color',[0.3 0.3 0.3]);
set(plot1(5),'DisplayName','V1','LineStyle','--','Color',[0.4 0.4 0.4]);
set(plot1(6),'DisplayName','V2','LineStyle',':','Color',[0.5 0.5 0.5]);
set(plot1(7),'DisplayName','V3','LineStyle','-.','Color',[0.6 0.6 0.6]);

% 取消以下行的注释以保留坐标区的 X 范围
xlim(axes1,[1.65415820358779 1.66673731880027]);
% 取消以下行的注释以保留坐标区的 Y 范围
ylim(axes1,[-0.000554837669973488 0.00088057056899609]);
box(axes1,'on');
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'ContextMenu','XTick',[0 1 2 3 4 5]);
% createaxes(axes2, t, Err_trace);