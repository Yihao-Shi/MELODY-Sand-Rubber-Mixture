%set up parameters
global x_min x_max y_min y_max;
global rigid_pos deformable_pos pos_vec;     
rigid_num = 1000;                                  %set required rigid balls number
deformable_num = 1000;                             %set required deformable balls number
tries_num = 200;                                 %set tries maximum number to generate balls
x_min = -2;
x_max = 2;
y_min = -1;
y_max = 0; 
rlo = 0.5;                                      %ball low radius
rhi = 0.5;                                      %ball high radius
rigid_pos = zeros(rigid_num,3);
deformable_pos = zeros(deformable_num,3);
particle_num = rigid_num+deformable_num;
pos_vec = zeros(particle_num,3);
overlap_flag = 0;
ball_count = 0;


%make assembly and check overlap
for i = 1:particle_num
    for j = 1:tries_num
        rc = rlo+rand(1,1)*(rhi-rlo);
        diameter=rc*2.0;
        xl = x_max-x_min-diameter;
        yl = y_max-y_min-diameter;
        xc = xl*(rand(1,1)+(x_min+rc)/xl);
        yc = yl*(rand(1,1)+(y_min+rc)/yl);
        
        %check overlap
        if ball_count>=1
            for m = 1:size(pos_vec,1)
                dist = sqrt(power((xc-pos_vec(m,1)),2)+power((yc-pos_vec(m,2)),2));
                dist_min = rc+pos_vec(m,3);
                if dist<dist_min
                    overlap_flag = 1;
                    break;
                else
                    overlap_flag = 0;
                end
            end
        end
        %store valid data
        if overlap_flag==1
            continue;
        else
            ball_count = ball_count+1;
            pos_vec(ball_count,:) = [xc,yc,rc];
            break
        end
    end
    if overlap_flag==1||j==tries_num
        disp('------Fewer balls ('+string(ball_count)+'balls) generated than specified');
        break;
    elseif ball_count==particle_num
        disp('--- '+string(ball_count)+' balls successfully generated');
    end
end


%plot general view of initial state
for n=1:size(pos_vec,1)                         
    X=pos_vec(n,1)+pos_vec(n,3)*cos((0:1:360)'/180*pi);
    Y=pos_vec(n,2)+pos_vec(n,3)*sin((0:1:360)'/180*pi);
    plot(X,Y);
    hold on; 
end
axis equal;


%distribute rigid balls and deformable balls
global id_number rigid_row deformable_row;
rigid_row = 1;
deformable_row = 1;
for idp=1:size(pos_vec,1)
    if floor(idp/2)<(idp/2)
       deformable_pos(deformable_row,:) = [pos_vec(idp,1),pos_vec(idp,2),pos_vec(idp,3)];
       deformable_row = deformable_row+1;
    else
	   rigid_pos(rigid_row,:) = [pos_vec(idp,1),pos_vec(idp,2),pos_vec(idp,3)];
       rigid_row = rigid_row+1;   
    end
end


%delete invalid data
pos_vec(all(pos_vec==0,2),:) = [];
rigid_pos(all(rigid_pos==0,2),:) = [];
deformable_pos(all(deformable_pos==0,2),:) = [];
rigid_row = size(rigid_pos,1);
deformable_row = size(deformable_pos,1);
id_number = size(pos_vec,1);

%find 