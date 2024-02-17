function [Err_dis] = ObjFunReprojErrcoutour(x, K,event_cur, P_3d,n)
%OBJFUNREPROJERRLINE 此处显示有关此函数的摘要


    RT(1:3,1:3) = rodrigues([x(1) x(2) x(3)]);
    RT(1:3,4) = [x(4); x(5); x(6)];
   
  
    Box_proj=K * RT* P_3d;
    p_proj= [Box_proj(1,:)./Box_proj(3,:); (Box_proj(2,:)./Box_proj(3,:));ones(1,size(Box_proj,2))];
%     biy=p_proj(1:2,:)';
%     k = boundary(biy);
%     contourPoints = biy(k,:);
% 
%      [idx, centroids] = kmeans(biy(k,:), num);
% 
% % 选择每个聚类的中心点作为代表性的边缘点
    contourPoints = p_proj(1:2,:)';
%plot(edge_points(:,1),edge_points(:,2));


% 设置最近邻点数目

%    n = 10;  % 调整 n 的值以适应你的需求
    
    % 初始化法向量矩阵
    normals = zeros(size(contourPoints));
    
    % 寻找每个点的最近邻点并计算法向量
    for i = 1:size(contourPoints, 1)
        point = contourPoints(i,:);
        
        % 计算每个点到所有点的欧氏距离
        distances = vecnorm(contourPoints - point, 2, 2);
        
        % 找到距离最近的 n 个点（除了自身,n=11）
        [~, indices] = mink(distances, 11);
        indices = indices(indices~=i); % 排除自身点
        
        % 提取最近的 n 个点的坐标
        nearestPoints = contourPoints(indices,:);
        
        % 根据最近的 n 个点计算法向量
        centroid = mean(nearestPoints);
        covMatrix = cov(nearestPoints);
        [V, D] = eig(covMatrix);
        [~, minIndex] = min(diag(D));
        normal = V(:, minIndex);
        
    %     % 判断法向量朝向，并反转方向使其朝向外侧
    %     if dot(normal, point - centroid) < 0
    %         normal = -normal;
    %     end
    %     
    %     % 统一法向量朝向外侧
    %     if dot(normal, point) < 0
    %         normal = -normal;
    %     end
        
        % 将法向量保存在法向量矩阵中
        normals(i,:) = normal';
    end


% 
%       hold on;
%      plot(contourPoints(:,1),contourPoints(:,2), 'r.');  % 显示平滑后的边缘点 
%        hold on
%          plot(event_cur(:,2),event_cur(:,3), 'b.');  % 显示平滑后的边缘点 %  
%      hold on
%      quiver(contourPoints(:,1),contourPoints(:,2),normals(:,1),normals(:,2), 'Color', 'g');  % 显示法向量
%             
%         hold on 
%            %  plot(corresponding_points(:,2),corresponding_points(:,1), 'black.'); 
%     hold off;


    [corresponding_points] = epcorres(contourPoints, normals, [event_cur(:,2),event_cur(:,3)],n);
    
   
      
    
    % 初始化存储距离的向量
 
    Err_dis = zeros(size(contourPoints,1), 1);
    
    % 对于每个边缘点和对应点
    for i = 1:size(contourPoints,1)
        % 检查对应点是否为空
        if ~(corresponding_points(i,1)==0 && corresponding_points(i,2)==0)
            % 进行投影并计算距离
            Err_dis(i) = norm(contourPoints(i,1:2) - corresponding_points(i,:));
        else
            Err_dis(i)= 0;
           
        end
    end

%  Err_dis(Err_dis == 0) = [];

end


function avgDistance = findAverageDistance(a, aNormals, b)
    n = size(a, 1);
    m = size(b, 1);
    distanceSum = 0;
    count = 0;
    
    for i = 1:n
        currentA = a(i,:);
        currentNormal = aNormals(i,:);
        
        angleRange = deg2rad(2);
        minAngle = atan2(currentNormal(2), currentNormal(1)) - angleRange;
        maxAngle = atan2(currentNormal(2), currentNormal(1)) + angleRange;
        
        % 计算扇形区域的单位向量
        fanVector = [cos(minAngle), sin(minAngle)];
        
        % 遍历每个点b
        for j = 1:m
            currentB = b(j,:);
            
            % 计算当前点a到点b的向量
            vectorAB = currentB - currentA;
            
            % 计算当前点a到点b的距离
            distance = norm(vectorAB);
            
            % 判断点b是否在扇形区域内
            if dot(vectorAB, fanVector) >= 0 && dot(vectorAB, currentNormal) >= 0 && acos(dot(vectorAB, currentNormal) / (norm(vectorAB) * norm(currentNormal))) <= angleRange
                distanceSum = distanceSum + distance;
                count = count + 1;
            end
        end

       avgDistance(i,1) = distanceSum / count;
    end
    
    
end
