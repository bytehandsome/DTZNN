function axes1 = createfigure_joint(X1, YMatrix1)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  x 数据的向量
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 11-May-2022 00:03:56 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','$$\theta_{1}$$', 'Color',[0 0 0]);
set(plot1(2),'DisplayName','$$\theta_{2}$$','LineStyle','--','Color',[0.1 0.1 0.1]);
set(plot1(3),'DisplayName','$$\theta_{3}$$','LineStyle',':','Color',[0.2 0.2 0.2]);
set(plot1(4),'DisplayName','$$\theta_{4}$$','LineStyle','-.','Color',[0.3 0.3 0.3]);
set(plot1(5),'DisplayName','$$\theta_{5}$$','LineStyle','--','Color',[0.4 0.4 0.4]);
set(plot1(6),'DisplayName','$$\theta_{6}$$','LineStyle',':','Color',[0.5 0.5 0.5]);
set(plot1(7),'DisplayName','$$\theta_{7}$$','LineStyle','-.','Color',[0.6 0.6 0.6]);

% 创建 ylabel
ylabel('关节角度/rad');

% 创建 xlabel
xlabel('时间/（秒）');

% 创建 title
% title('joint__trace');

box(axes1,'on');
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'XTick',[0 1 2 3 4 5]); % ,'YTick',[-3 -2 -1 0 1 2 3]
% 创建 legend
legend1 = legend(axes1,'show','Location','Best');
set(legend1,'Interpreter','latex');
% axes1 = createfigure_joint(t, q_trace);
% saveas(gcf,'joint','png');