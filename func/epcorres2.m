%矩形区域搜索
function [corresponding_points] = epcorres2(edge_points, normals, n_points,n)
    corresponding_points = [];

    for i = 1:size(edge_points, 1)
        edge_point = edge_points(i, :);
        normal = normals(i, :)/norm(normals(i, :));







vec = bsxfun(@minus, n_points, edge_point); % 计算每个点到直线上点的向量
distances1 = abs(vec * normal' / norm(normal)); % 计算每个点到直线的距离

vec2 = bsxfun(@minus, n_points, edge_point); % 计算每个点到直线上点的向量
distances2 = sqrt(sum(vec .^ 2, 2)); % 计算每个点到直线的距离

indices = find(distances1 <= n & distances2 <= n*sqrt(2));
points_within_threshold = n_points(indices, :);
        
 
        corresponding_points_found = points_within_threshold;
        % 计算距离在10个像素以内的所有点的中位数的点
        if ~isempty(corresponding_points_found)
          corresponding_point = mean(corresponding_points_found,1);
       %       corresponding_point=kmeans(corresponding_points_found, 1);
        else
            corresponding_point = [0,0]; % 如果没有找到点，则设置为 [0,0]
        end
        
        corresponding_points = [corresponding_points;corresponding_point];
    end


%       plotCorrespondingPoints(edge_points, corresponding_points);
 
  
  
%      quiver(edge_points(i, 2), edge_points(i, 1),normal(2), normal(1), 'Color', 'b');  % 显示法向量
%        hold on
%     scatter(corresponding_points_found(:, 2), corresponding_points_found(:, 1), 'red');

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