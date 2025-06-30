%set up parameters
global x_min x_min_expand x_max x_max_expand y_min y_min_expand y_max y_max_expand rhi rlo particle_num;
global rigid_pos deformable_pos pos_vec;   
pos_vec=load('./ball_pos_PFC/Ball_pos.txt','r');
info=textread('./ball_pos_PFC/Ball&Wall_info.txt','%f');
x_min=info(1);
x_max=info(2);
y_min=info(3);
y_max=info(4);
x_min_expand = 1.5*x_min;
x_max_expand = 1.5*x_max;
y_min_expand = 1.5*y_min;
y_max_expand = 1.5*y_max;
rlo=info(5);
rhi=info(6);
particle_num=info(7);
proportion = 0;                             %the proportion of soft grains in mixture
rigid_num = ceil((1-proportion)*particle_num);
deformable_num = particle_num-rigid_num;
rigid_pos = zeros(rigid_num,3);
rigid_row = 1;
deformable_row = 1;

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