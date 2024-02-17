%圆形区域搜索
function [corresponding_points] = epcorres(edge_points, normals, scattered_points,n)
    corresponding_points = zeros(size(edge_points, 1), 2);
    
    for i = 1:size(edge_points, 1)
        edge_point = edge_points(i, :);
        normal = normals(i, :);

        distances = sqrt(sum((scattered_points - edge_point).^2, 2));

% 找到距离小于5的所有点
        corresponding_points_found = scattered_points(distances < n, :);
%         % 根据法向量计算对应点坐标
%         corresponding_point_pos = edge_point + normal;
%         corresponding_point_neg = edge_point - normal;
%         
%         % 沿着正向法向量方向搜索，并找到距离不超过10个像素的所有点
%         pos_indices = vecnorm(scattered_points - edge_point, 2, 2) <= n;
%         corresponding_points_found_pos = scattered_points(pos_indices, :);
% 
%         corresponding_points_found = [corresponding_points_found_pos];

        % 计算距离在10个像素以内的所有点的中位数的点
        if ~isempty(corresponding_points_found)
            corresponding_point = mean(corresponding_points_found, 1);
           %[Inx  corresponding_point]=kmeans(corresponding_points_found, 1);
        else
            corresponding_point = [0,0]; % 如果没有找到点，则设置为 [0,0]
        end



        
        corresponding_points(i,:) = corresponding_point;
    end


%     scatter(corresponding_points_found_pos (:,1),corresponding_points_found_pos (:,2),'red')
%     hold on
%     scatter(corresponding_points_found_neg (:,1),corresponding_points_found_neg (:,2),'blue')
%     hold on
%      quiver(edge_point(:,1),edge_point(:,2),normal(:,1),normal(:,2), 'Color', 'g');  % 显示法向量
%      hold on
%      scatter(corresponding_points (i,1),corresponding_points (i,2),'black')
%      hold off;

    % 可以取消下面这行的注释，用于绘制边缘点和对应点
%    plotCorrespondingPoints(edge_points, corresponding_points)
%     plotqumian(corresponding_points_found)
end

function plotCorrespondingPoints(edge_points, corresponding_points)
    figure;
    
    % 绘制边缘点
    scatter(edge_points(:, 1), edge_points(:, 2), 'b');
    hold on;
    
    % 绘制对应点
    for i = 1:size(corresponding_points, 1)
        point = corresponding_points(i,:);
        if ~isempty(point)
            scatter(point(1), point(2), 'r');
        end
    end
    
    % 设置图例和标题
    legend('Edge Points', 'Corresponding Points');
    title('Corresponding Points');
    
    hold off;
end


function plotqumian(corresponding_points_found)
[unique_points, ~, idx] = unique(corresponding_points_found, 'rows');
counts = accumarray(idx, 1);

% 更新 corresponding_points_found 的第三列为出现次数
corresponding_points_found(:, 3) = counts(idx);

corresponding_points = unique(corresponding_points_found, 'rows');
x = corresponding_points(:,1);
y =corresponding_points(:,2);
z = corresponding_points(:,3);
% % 生成网格
% [X,Y] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
% 
% % 进行插值
% Z = griddata(x, y, z, X, Y, 'cubic');
% 
% % 绘制曲面
% surf(X, Y, Z);
% 
% % 绘制三维曲面图
% surf(Xq, Yq, Vq, 'EdgeColor', 'none');
% colormap(jet);
% 
% % 设置坐标轴标签和图标题
% xlabel('X Label');
% ylabel('Y Label');
% zlabel('Z Label');
% title('Corresponding Points Found Density');
% 
% % 标注最高点
% hold on;
% plot3(max_x, max_y, max_z, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% text(max_x, max_y, max_z, ['Max: (' num2str(max_x) ',' num2str(max_y) ',' num2str(max_z) ')'], 'FontSize', 12);
% hold off;
% 统计每个坐标点出现的次数
% x = corresponding_points_found(:, 1);
% y = corresponding_points_found(:, 2);
% [unique_points, ~, idx] = unique([x y], 'rows');
% counts = accumarray(idx, 1);

% 绘制三维柱状图
bar3([x,y,z],0.7);
xlabel('X');
ylabel('Y');
zlabel('Counts');
title('Corresponding Points');

% 设置图形美化选项

end