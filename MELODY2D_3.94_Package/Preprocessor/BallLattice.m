%set up parameters
global x_min x_max y_min y_max rhi rlo particle_num;
global rigid_pos deformable_pos pos_vec;     
row_num = 3;                                %set required rigid balls number
column_num = 3;                             %set required deformable balls number 
rhi = 0.01;                                 %the maximum of ball radius
rlo = 0.01;                                 %the minimum of ball radius
x_min = -(column_num+0.5)*rhi;
x_max = (column_num+0.5)*rhi;
y_min = -0.5*((row_num-1)*sqrt(3)+2)*rhi;
y_max = 0.5*((row_num-1)*sqrt(3)+2)*rhi;
particle_num = row_num*column_num;
proportion = 0;                             %the proportion of soft grains in mixture
rigid_num = ceil((1-proportion)*particle_num);
deformable_num = particle_num-rigid_num;
rigid_pos = zeros(rigid_num,3);
deformable_pos = zeros(deformable_num,3);
pos_vec = zeros(particle_num,3);
ball_count = 0;
rigid_row = 1;
deformable_row = 1;

%place the particles
for i = 1:row_num
    for j = 1:column_num
        ball_count = ball_count+1;
        if floor(i/2)<i/2
            pos_vec(ball_count,1) = x_min+2*(rhi+(j-1)*rhi);
        else
            pos_vec(ball_count,1) = x_min+rhi+2*(j-1)*rhi;
        end
        pos_vec(ball_count,2) = y_min+rhi+sqrt(3)*(i-1)*rhi;
        pos_vec(ball_count,3) = rlo+rand(1,1)*(rhi-rlo);            %ball radius;
    end
end

%Update x_length of shear box
x_min = -6*ceil((column_num+0.5)/6)*rhi;
x_max = 6*ceil((column_num+0.5)/6)*rhi;

%plot general view of initial state
for n=1:size(pos_vec,1)                         
    X=pos_vec(n,1)+pos_vec(n,3)*cos((0:1:360)'/180*pi);
    Y=pos_vec(n,2)+pos_vec(n,3)*sin((0:1:360)'/180*pi);
    plot(X,Y);
    hold on; 
end
axis equal;

%distribute rigid balls and deformable balls
pos_copy = pos_vec;
for c = 1:rigid_num
    pos_rand = ceil(size(pos_copy,1)*rand(1,1));
    rigid_pos(rigid_row,:) = pos_copy(pos_rand,:);
    rigid_row = rigid_row+1;
    pos_copy(pos_rand,:) = [];
end
deformable_pos = pos_copy;

%output ball information
% xlswrite('Position_info.xlsx',[x_min,x_max;y_min,y_max],'Overall_info');
% xlswrite('Position_info.xlsx',pos_vec,'Overall_info','C1');
% xlswrite('Position_info.xlsx',rigid_pos,'Rigid_info');
%xlswrite('Position_info.xlsx',deformable_pos,'Deformable_info');