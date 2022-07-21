function  createfigure_err(X1, YMatrix1)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  x 数据的向量
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 11-May-2022 10:29:44 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','$$e_{x}$$ ( m )', 'Color',[0 0 0]);
set(plot1(2),'DisplayName','$$e_{y}$$ ( m )','LineStyle','--','Color',[0.1 0.1 0.1]);
set(plot1(3),'DisplayName','$$e_{z}$$ ( m )','LineStyle',':','Color',[0.2 0.2 0.2]);
set(plot1(4),'DisplayName','1-$$q_{s}$$(rad)','LineStyle','-.','Color',[0.3 0.3 0.3]);
set(plot1(5),'DisplayName','$$q_{v1}$$ (rad)','LineStyle','--','Color',[0.4 0.4 0.4]);
set(plot1(6),'DisplayName','$$q_{v2}$$ (rad)','LineStyle',':','Color',[0.5 0.5 0.5]);
set(plot1(7),'DisplayName','$$q_{v3}$$ (rad)','LineStyle','-.','Color',[0.6 0.6 0.6]);

% 创建 ylabel
ylabel('误差');

% 创建 xlabel
xlabel('时间/（秒）');

% 创建 title
% title('Err__trace');

box(axes1,'on');
hold(axes1,'on');
% 设置其余坐标区属性
set(axes1,'XTick',[0 1 2 3 4 5]);
% 创建 legend
legend1 = legend(axes1,'show','Location','Best');
set(legend1,'Interpreter','latex');

% 创建 axes
axes2 = axes('Parent',figure1,...
    'Position',[0.317222222222221 0.255358565737051 0.309999999999999 0.31]);
hold(axes2,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1,'Parent',axes2,'LineWidth',2);
set(plot1(1),'DisplayName','$$e_{x}$$', 'Color',[0 0 0]);
set(plot1(2),'DisplayName','$$e_{y}$$','LineStyle','--','Color',[0.1 0.1 0.1]);
set(plot1(3),'DisplayName','$$e_{z}$$','LineStyle',':','Color',[0.2 0.2 0.2]);
set(plot1(4),'DisplayName','$$1-q_{s}$$','LineStyle','-.','Color',[0.3 0.3 0.3]);
set(plot1(5),'DisplayName','$$q_{v1}$$','LineStyle','--','Color',[0.4 0.4 0.4]);
set(plot1(6),'DisplayName','$$q_{v2}$$','LineStyle',':','Color',[0.5 0.5 0.5]);
set(plot1(7),'DisplayName','$$q_{v3}$$','LineStyle','-.','Color',[0.6 0.6 0.6]);

% 取消以下行的注释以保留坐标区的 X 范围
xlim(axes2,[3.8  4.2]);
% 取消以下行的注释以保留坐标区的 Y 范围
ylim(axes2,[-0.00015 0.00015]);
box(axes2,'on');
hold(axes1,'off');
% 设置其余坐标区属性
% set(axes2,'ContextMenu','XTick',[0 1 2 3 4 5]);
% set(axes2,'ContextMenu','YTick',[-1 0 1]);

% createfigure_err(t, Err_trace);
% saveas(gcf,'err','png');