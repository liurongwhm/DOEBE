function gbest = FindGlobalBest(GBA,f_GBA,x,f_x)
% Find gbest from the Global Best Archive£¨GBA£©for each solution in x
[x_num,dim] = size(x);
[GBA_num,~] = size(GBA);
gbest = zeros(x_num,dim);
k = max(3,round(GBA_num / x_num));

temp1 = f_x(:,1);
temp2 = f_x(:,2);
x_sigma = (temp1.^2 - temp2.^2) ./ (temp1.^2 + temp2.^2);

temp1 = f_GBA(:,1);
temp2 = f_GBA(:,2);
GBA_sigma = (temp1.^2 - temp2.^2) ./ (temp1.^2 + temp2.^2);

temp = abs(ones(GBA_num,1)*x_sigma'-GBA_sigma*ones(1,x_num));
[val,idx] = sort(temp); %sort each column

CrowdingDistance = zeros(GBA_num,1);
CrowdingDistance(1) =  1000000;
CrowdingDistance(GBA_num) = 1000000;
m1 =  f_GBA(GBA_num,1) - f_GBA(1,1);
m2 = f_GBA(1,2) - f_GBA(GBA_num,2);

for i = 2:GBA_num-1
    CrowdingDistance(i) = (f_GBA(i+1,1) -  f_GBA(i-1,1))/m1 + (f_GBA(i-1,2) -  f_GBA(i+1,2))/m2;
end

if GBA_num<5
    for i =1:x_num
        gbest(i,:) = GBA(randperm(GBA_num,1),:);
    end
else
    idx = idx(1:k,:);
    for i = 1:x_num
        index = idx(:,i);
        temp3 = CrowdingDistance(index,1);
        [~,temp4] = max(temp3);
        gbest(i,:) = GBA(index(temp4),:);
    end
end



end

