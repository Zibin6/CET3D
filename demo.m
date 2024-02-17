
%尝试进行event的位姿估计 2023/12/18 success
close all
  clc
  clear
addpath('func\')
     %fileD=fopen('.\result_crs.txt','w+');

load("data.mat")
  

image_num=1;


K=[320	0	320
0	320	240
0	0	1 ];



% % 去畸变的代码
event_t=0.02;
enum=2000;
event_number=5000;
% for k=1:event_number:size(event,1)  % 主循环
for k=event(1,1):event_t:event(end,1)
 
%     p_cur=K*pose_cur*P_3d;
%     p_cur= [p_cur(1,:)./p_cur(3,:); (p_cur(2,:)./p_cur(3,:));ones(1,size(p_cur,2))];
    
% clustering 
    [~, index] = min(abs(k - event(:,1)));
    start_index = index(1) - enum;
    end_index = index(1) + enum;
    start_index=max(1,start_index);

    event_cur=event(start_index:end_index,:);

    event_proj=K *pose_cur* object;
    p_proj= [event_proj(1,:)./event_proj(3,:); (event_proj(2,:)./event_proj(3,:));ones(1,size(event_proj,2))];
    biy=p_proj(1:2,:)';
    aa = boundary(biy);

    object1=object(:,aa)';
   

    [idx, contour] = kmeans(object1, 200);

    x_init = zeros(6, 1);
    ang_init = rodrigues(pose_cur(1:3,1:3));
    x_init(1:3) = ang_init;
    x_init(4:6) = pose_cur(1:3,4);


    options = optimset;
    options = optimoptions("lsqnonlin","Display","none");
    options.Algorithm = 'levenberg-marquardt';
    options.MaxFunEvals = 20000;
    options.TolFun = 1e-5;
    options.TolX = 1e-5;
    options.MaxIter = 1000;

    [x_optim,resnorm,residual] =lsqnonlin(@(x) ObjFunReprojErrcoutour( x, K,event_cur, contour',5), x_init, [], [], options);
    
    R_opt = rodrigues([x_optim(1) x_optim(2) x_optim(3)]);

    T_opt = [x_optim(4); x_optim(5); x_optim(6)];

    pose_cur=[R_opt T_opt];
   
    q=dcm2quat(pose_cur(1:3,1:3));

   % fprintf(fileD, '%.6f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n', k,T_opt(1), T_opt(2), T_opt(3),q(2), q(3), q(4), q(1) );

    repro_Image=reproj_show122(K,pose_cur,event_cur,object1');


end


function [repro_Image]=reproj_show122(K,pose_cur,event_cur,P_3d)

% P_3d 4*n,第四行是1


 p_pro=K*[pose_cur]*P_3d;
 p_pro= [p_pro(1,:)./p_pro(3,:); p_pro(2,:)./p_pro(3,:);ones(1,size(p_pro,2))];
 image_cur=zeros(480,640);
    for i=1:size(event_cur,1)
        if event_cur(i,3)==0 
            event_cur(i,3)=1;
        end
        if event_cur(i,2)==0 
            event_cur(i,2)=1;
        end
      image_cur(event_cur(i,3), event_cur(i,2)) = 255;
    end
    imshow(image_cur);


     hold on
     plot(p_pro(1,:),p_pro(2,:),'color','r','LineWidth',1);
     line([pp1(1,:);pp2(1,:)],[pp1(2,:);pp2(2,:)],'color','r','LineWidth',1);
     hold off
   
    repro_Image = getframe;
    repro_Image = repro_Image.cdata;


end

