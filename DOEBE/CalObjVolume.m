function [objs] = CalObjVolume(IndexEndmember,X)
% volume calculation

[PopNum,p] = size(IndexEndmember);
SubNum = size(X,2);
objs = zeros(PopNum,SubNum);

%Calculate the objective function values
for i = 1:PopNum
    for j = 1:SubNum
        % Calculate the inverse volume of the simplex constructed by endmembers
        Xtemp = X{j};
        EndmemberReduction = Xtemp(:,IndexEndmember(i,:));
        objs(i,j) = 1/abs(det(EndmemberReduction));      
    end
%     EndmemberReduction = imgTrans(:,IndexEndmember(i,:));
%     objs(i,j+1) = 1/abs(det(EndmemberReduction));      
end
    
end

