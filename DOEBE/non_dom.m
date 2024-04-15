function [newGBA,newf_GBA] = non_dom(oldGBA,oldf_GBA,pbest,f_pbest)
% Input£º f_x£¬Each column represents the objective values of all particles
% Output£ºSorted non dominant solutions
newGBA = [];
newf_GBA = [];
x = [oldGBA;pbest];
f_x = [oldf_GBA;f_pbest];
[r,c] = size(f_x);
f1 = f_x(:,1);
f2 = f_x(:,2);

for i = 1:r
    temp = ((f1<f1(i))+(f2<f2(i)));
    if max(temp) <2
        newGBA = [newGBA;x(i,:)];
        newf_GBA = [newf_GBA;f_x(i,:)];
    end
end

[newGBA,ia,~] = unique(newGBA,'rows');
newf_GBA = newf_GBA(ia,:);

[newf_GBA(:,1),idx] = sort(newf_GBA(:,1));
newf_GBA(:,2) = newf_GBA(idx,2);
newGBA = newGBA(idx,:);


end

