function axes1 = createfigure_qdot(X1, YMatrix1)
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
set(plot1(1),'DisplayName','$$\dot{\theta}_{1}$$', 'Color',[0 0 0]);
set(plot1(2),'DisplayName','$$\dot{\theta}_{2}$$','LineStyle','--','Color',[0.1 0.1 0.1]);
set(plot1(3),'DisplayName','$$\dot{\theta}_{3}$$','LineStyle',':','Color',[0.2 0.2 0.2]);
set(plot1(4),'DisplayName','$$\dot{\theta}_{4}$$','LineStyle','-.','Color',[0.3 0.3 0.3]);
set(plot1(5),'DisplayName','$$\dot{\theta}_{5}$$','LineStyle','--','Color',[0.4 0.4 0.4]);
set(plot1(6),'DisplayName','$$\dot{\theta}_{6}$$','LineStyle',':','Color',[0.5 0.5 0.5]);
set(plot1(7),'DisplayName','$$\dot{\theta}_{7}$$','LineStyle','-.','Color',[0.6 0.6 0.6]);


% 创建 ylabel
ylabel('关节速度/（弧度/秒）');

% 创建 xlabel
xlabel('时间/（秒）');

% 创建 title
% title('joint__trace');

box(axes1,'on');
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'XTick',[0 1 2 3 4 5],'YTick',[-0.5 0 0.5 1 1.5 2]);
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,'Interpreter','latex');

% 创建 axes
axes2 = axes('Parent',figure1,...
    'Position',[0.332698412698411 0.401390311768797 0.309999999999999 0.31]);
hold(axes2,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1,'Parent',axes2,'LineWidth',2);
set(plot1(1),'DisplayName','$$\dot{\theta}_{1}$$', 'Color',[0 0 0]);
set(plot1(2),'DisplayName','$$\dot{\theta}_{2}$$','LineStyle','--','Color',[0.1 0.1 0.1]);
set(plot1(3),'DisplayName','$$\dot{\theta}_{3}$$','LineStyle',':','Color',[0.2 0.2 0.2]);
set(plot1(4),'DisplayName','$$\dot{\theta}_{4}$$','LineStyle','-.','Color',[0.3 0.3 0.3]);
set(plot1(5),'DisplayName','$$\dot{\theta}_{5}$$','LineStyle','--','Color',[0.4 0.4 0.4]);
set(plot1(6),'DisplayName','$$\dot{\theta}_{6}$$','LineStyle',':','Color',[0.5 0.5 0.5]);
set(plot1(7),'DisplayName','$$\dot{\theta}_{7}$$','LineStyle','-.','Color',[0.6 0.6 0.6]);

% 取消以下行的注释以保留坐标区的 X 范围
xlim(axes2,[3.8  4.2]);
% 取消以下行的注释以保留坐标区的 Y 范围
% ylim(axes2,[-0.000254837669973488 0.00028057056899609]);
box(axes2,'on');
% hold(axes1,'off');
% 设置其余坐标区属性
% set(axes2,'ContextMenu','XTick',[0 1 2 3 4 5]);



% createfigure_qdot(t, jointVelocity);
% saveas(gcf,'qdot','jpg');

