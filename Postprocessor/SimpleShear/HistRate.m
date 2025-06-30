function [result]=HistRate(force)
    result=zeros(22,2);
    result(:,1)=[0,10.^((-300:20:100)/100)].';
    force=sort(force);
    for i=2:size(result,1)
        for j=result(i-1,2)+1:length(force)
            if force(j)>=result(i-1,1)&force(j)<result(i,1)
                result(i,2)=result(i,2)+1;
            end
        end
    end
    result(:,2)=result(:,2)/length(force);
end