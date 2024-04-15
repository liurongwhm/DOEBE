function [indicies,record] = DOEBE(img2d,img3d,SubNum,PopNum1,PopNum2,P,r,minPixels)
% A Novel Endmember Bundle Extraction Framework for Capturing Endmember Variability by Dynamic Optimization（DOEBE）
% Input:
%   img3d:      3-D matrix of the image data (row×col×L)
%   img2d:      2-D matrix of the image data (L×N)
%   P:          Number of endmembers
%   Subnum:     Number of subimage (Default is 2)
%   PopNum1:    Population size of the multi-objective optimization (Default is 20)
%   PopNum2:    Population size of the single-objective optimization (Default is 20)
%   r:          Threshold of the reconstruction error (Default is 0.006)
%   minpixel:   Minimum number of pixels required for the next iteration (Default is 0.01*N)
% Output:
%   indicies:   Extracted endmembers
%   record:     Minimum objective value and reconstructed error of remaining pixels of each iteration

[row,col,L] = size(img3d);
[~,N] = size(img2d);
img2d_reduced = img2d;
Pixels = N;
proM = 1;  % Mutation rate of the polynomial mutation
disM = 20; % Distribution index of the polynomial mutation
indicies = [];
record = [];

%% Step1：Population initialization based on MODPSO
subImg2d = cell(1,SubNum);
X = cell(1,SubNum);

% Reduce the dimension of the hyperspectral image
[~,A] = hyperMnf(img2d, row, col);
A = A';
trans = A(1:P-1,:);
SubTrans = cell(1,SubNum);
imgTrans = [ones(1,N);trans * img2d];

% Equally divide the original image into two subsets along the spectral dimension
SubBands = floor(L/SubNum);
for i = 1:SubNum-1
    subImg2d{1,i} = img2d((i-1)*SubBands+1:i*SubBands,:);
    SubTrans{1,i} = trans(:,(i-1)*SubBands+1:i*SubBands);
    X{1,i} = [ones(1,N);SubTrans{1,i} * subImg2d{1,i}];
end
subImg2d{1,SubNum} = img2d((SubNum-1)*SubBands+1:end,:);
SubTrans{1,SubNum} = trans(:,(SubNum-1)*SubBands+1:end);
X{1,SubNum} = [ones(1,N);SubTrans{1,SubNum} * subImg2d{1,SubNum}];

% Initialize the population
x = zeros(PopNum1,N);  % solutions
for j=1:PopNum1
    x(j,randperm(N,P)) = 1;
    IndexEndmember(j,:) = find(x(j,:)==1);
end
pbest = x;
f_x = CalObjVolume(IndexEndmember,X);
f_pbest = f_x;
[GBA,f_GBA] = non_dom([],[],pbest,f_pbest);

% Multiobjective discrete particle swarm optimization (MODPSO)
maxiter = 300;
iter_finish = 0;
prob = 0.2;  % Probability of random movement of particles
iter = 0;
GBA_record1= zeros(maxiter,1);
GBA_record2= zeros(maxiter,1);
while iter < maxiter && iter_finish == 0
    iter = iter+1;
    gbest = FindGlobalBest(GBA,f_GBA,x,f_x);
    % Update speed and current position
    for j=1:PopNum1
        % Update speed
        if rand()<prob % Move randomly
            zero_idx = find(x(j,:)==0);
            one_idx = find(x(j,:)==1);
            v_positive_idx = zero_idx(randperm(length(zero_idx),1));
            v_negative_idx = one_idx(randperm(length(one_idx),1));
        else % Move according to pbest and gbest
            vp = pbest(j,:) - x(j,:);
            vg = gbest(j,:) - x(j,:);
            vs= vp+vg;
            positive_idx = find(vs>0);
            negative_idx = find(vs<0);
            if ~isempty(positive_idx)
                v_positive_idx = positive_idx(randperm(length(positive_idx),1));
                v_negative_idx = negative_idx(randperm(length(negative_idx),1));
            else
                zero_idx = find(x(j,:)==0);
                one_idx = find(x(j,:)==1);
                v_positive_idx = zero_idx(randperm(length(zero_idx),1));
                v_negative_idx = one_idx(randperm(length(one_idx),1));
            end
        end
        % Update position
        x(j,v_positive_idx)=1;
        x(j,v_negative_idx)=0;
        IndexEndmember(j,:) = find(x(j,:)==1);
    end
    % Update pbest
    f_x = CalObjVolume(IndexEndmember,X);
    for j=1:PopNum1
        n = ((f_x(j,1)<f_pbest(j,1)) + (f_x(j,2)<f_pbest(j,2)));
        if n==2
            pbest(j,:) = x(j,:);
            f_pbest(j,:) = f_x(j,:);
        end
        if n==1
            if rand() > 0.5
                pbest(j,:) = x(j,:);
                f_pbest(j,:) = f_x(j,:);
            end
        end
    end
    [GBA,f_GBA] = non_dom(GBA,f_GBA,pbest,f_pbest);
    f_GBA
    GBA_record1(iter)= min(f_GBA(:,1));
    GBA_record2(iter)= min(f_GBA(:,2));
    fprintf('Excecuting MODPSO，the %dth iteration finished\n',iter);
end
record{1,1} = GBA_record1;
record{1,2} = GBA_record2;

%% Step2：Endmember extraction based on CCSO
k = 0;
kmax = 30;
while Pixels >= minPixels && k < kmax
    k = k+1;
    if r<0.02
        r = r+(k-1)*0.002;
    end
    % Initialization
    x1 = zeros(PopNum2,2*P);
    x2 = x1;
    for j = 1:PopNum2
        x1(j,:) = [randperm(row,P) randperm(col,P)];
        x2(j,:) = [randperm(row,P) randperm(col,P)];
    end
    % Initialize the population with the results of MODPSO
    x1_binary(1:round(PopNum2/2),:) = x(randperm(PopNum1,round(PopNum2/2)),:);
    x2_binary(1:round(PopNum2/2),:) = x(randperm(PopNum1,round(PopNum2/2)),:);
    for i = 1:round(PopNum2/2)
        x1_real = find(x1_binary(i,:)==1);
        x2_real = find(x2_binary(i,:)==1);
        for j = 1:P
            C = floor((x1_real(j)-1)/row)+1;
            R = x1_real(j)-row*(C-1);
            x1(i,j+P) = C;
            x1(i,j) = R;
            C = floor((x2_real(j)-1)/row)+1;
            R = x2_real(j)-row*(C-1);
            x2(i,j+P) = C;
            x2(i,j) = R;
        end
    end
    for j=1:PopNum2
        binarycode = transformRCToBinary(row,col,P,x1(j,:));
        f_x1(j,1) = unmixed(img2d_reduced,img2d(:,binarycode==1),1);
        binarycode = transformRCToBinary(row,col,P,x2(j,:));
        f_x2(j,1) = unmixed(img2d_reduced,img2d(:,binarycode==1),1);
    end
    xp = x1;
    f_xp = f_x1;
    v = zeros(PopNum2,2*P);
    v2 = [];
    iter = 0;
    iter_finish = 0;
    % CCSO
    while iter < maxiter && iter_finish == 0
        iter = iter+1;
        %% Competitive Update
        xx = [];
        tempx1 = x1;
        tempv = v;
        while size(tempx1,1) > 1
            select_id = randperm(size(tempx1,1),2);
            p = tempx1(select_id(1),:);
            q = tempx1(select_id(2),:);
            tempx1(select_id,:) = [];
            if f_x1(select_id(1)) > f_x1(select_id(2))
                xw = q;
                xl = p;
                wid = select_id(2);
                lid = select_id(1);
            else
                xw = p;
                xl = q;
                wid = select_id(1);
                lid = select_id(2);
            end
            r0 = rand();
            r1 = rand();
            sigma = 2*2*(1-iter/maxiter).^2+1;
            v2 = [v2;r0*rand()*tempv(lid,:)+sigma*r1*(xw-xl);tempv(wid,:)];
            xl = xl+v2(end-1,:)+r0*(v2(end-1,:)-tempv(lid,:));
            tempv(select_id,:) = [];
            off = [xl;xw];
            % Polynomial mutation
            n = 1;
            lower = ones(1,2*P);
            upper = [row*ones(1,P) col*ones(1,P)];
            Lower = repmat(lower,2*n,1);
            Upper = repmat(upper,2*n,1);
            Site  = rand(2*n,2*P) < proM/(2*P);
            mu    = rand(2*n,2*P);
            temp  = Site & mu<=0.5;
            off       = min(max(off,Lower),Upper);
            off(temp) = round(off(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                (1-(off(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1));
            temp = Site & mu>0.5;
            off(temp) = round(off(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                (1-(Upper(temp)-off(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))));
            off = round(off);
            % De-duplication
            for i = 1:size(off,1)
                binarycode = transformRCToBinary(row,col,P,off(i,:));
                positive_num = sum(binarycode);
                if positive_num < 6
                    temp = find(binarycode==1);
                    for j = 1:positive_num
                        C = floor((temp(j)-1)/row)+1;
                        R = temp(j)-row*(C-1);
                        off(i,j+P) = C;
                        off(i,j) = R;
                    end
                    for j = positive_num+1:P
                        off(i,j) = randperm(row,1);
                        off(i,j+P) = randperm(col,1);
                    end
                end
            end
            xx = [xx;off];
        end
        %% Selection
        for j = 1:PopNum2
            binarycode = transformRCToBinary(row,col,P,xx(j,:));
            f_xx(j,1) = unmixed(img2d_reduced,img2d(:,binarycode==1),1);
        end
        x1 = [x1;xx];
        f_x1 = [f_x1;f_xx];
        v = [v;v2];

        % De-duplication
        [x1,ia,~] = unique(x1,"rows");
        f_x1 = f_x1(ia,:);
        v = v(ia,:);

        [f_x1,I] = sort(f_x1);
        x1 = x1(I,:);
        v = v(I,:);
        if size(x1,1) > PopNum2
            x1 = x1(1:PopNum2,:);
            f_x1 = f_x1(1:PopNum2,:);
            v = v(1:PopNum2,:);
        end

        %% Cooperative Update
        % Learning pool selection
        Lp = [];
        while size(Lp,1) < PopNum2
            select_id = randperm(size(x2,1),2);
            p = x2(select_id(1),:);
            q = x2(select_id(2),:);
            if f_x2(select_id(1)) < f_x2(select_id(2))
                Lp = [Lp;p];
            elseif f_x2(select_id(1)) > f_x2(select_id(2))
                Lp = [Lp;q];
            else
                temp = [p;q];
                Lp = [Lp;temp(randperm(2,1),:)];
            end
        end
        y = [];
        while size(y,1) < PopNum2
            select_id = randperm(size(Lp,1),2);
            xw = Lp(select_id(1),:);
            xl = Lp(select_id(2),:);
            a = round(rand(1,2*P));
            xw = a.*xw+(1-a).*xl;
            xl = a.*xl+(1-a).*xw;
            off = [xl;xw];
            % Polynomial mutation
            n = 1;
            lower = ones(1,2*P);
            upper = [row*ones(1,P) col*ones(1,P)];
            Lower = repmat(lower,2*n,1);
            Upper = repmat(upper,2*n,1);
            Site  = rand(2*n,2*P) < proM/(2*P);
            mu    = rand(2*n,2*P);
            temp  = Site & mu<=0.5;
            off       = min(max(off,Lower),Upper);
            off(temp) = round(off(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                (1-(off(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1));
            temp = Site & mu>0.5;
            off(temp) = round(off(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                (1-(Upper(temp)-off(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))));
            off = round(off);
            % De-duplication
            for i = 1:size(off,1)
                binarycode = transformRCToBinary(row,col,P,off(i,:));
                positive_num = sum(binarycode);
                if positive_num < 6
                    temp = find(binarycode==1);
                    for j = 1:positive_num
                        C = floor((temp(j)-1)/row)+1;
                        R = temp(j)-row*(C-1);
                        off(i,j+P) = C;
                        off(i,j) = R;
                    end
                    for j = positive_num+1:P
                        off(i,j) = randperm(row,1);
                        off(i,j+P) = randperm(col,1);
                    end
                end
            end
            y = [y;off];
        end
        %% Selection
        for j = 1:PopNum2
            binarycode = transformRCToBinary(row,col,P,y(j,:));
            f_y(j,1) = unmixed(img2d_reduced,img2d(:,binarycode==1),1);
        end
        x2 = [x2;y];
        f_x2 = [f_x2;f_y];

        % De-duplication
        [x2,ia,~] = unique(x2,"rows");
        f_x2 = f_x2(ia,:);

        [f_x2,I] = sort(f_x2);
        x2 = x2(I,:);
        if size(x2,1) > PopNum2
            x2 = x2(1:PopNum2,:);
            f_x2 = f_x2(1:PopNum2,:);
        end
        %% Selection
        xp = [xp;xx;y];
        f_xp = [f_xp;f_xx;f_y];

        % De-duplication
        [xp,ia,~] = unique(xp,"rows");
        f_xp = f_xp(ia,:);

        [f_xp,I] = sort(f_xp);
        xp = xp(I,:);
        if size(xp,1) > 2*PopNum2
            xp = xp(1:2*PopNum2,:);
            f_xp = f_xp(1:2*PopNum2,:);
        end

        %% Inter-study
        f_gbest1 = min(f_x1);
        f_gbest2 = min(f_x2);
        gbest1 = x1(find(f_x1==f_gbest1),:);
        gbest2 = x2(find(f_x2==f_gbest2),:);
        if size(gbest1,1) > 1
            gbest1 = gbest1(randperm(size(gbest1,1),1),:);
        end
        if size(gbest2,1) > 1
            gbest2 = gbest2(randperm(size(gbest2,1),1),:);
        end
        if gbest1 < gbest2
            x2(end,:) = gbest1;
            f_x2(end,:) = f_gbest1;
        elseif gbest1 > gbest2
            x1(end,:) = gbest2;
            f_x1(end,:) = f_gbest2;
            v(end,:) = 0;
        end

        %% Find gbest
        f_gbest_xp = f_xp(1);
        gbest_xp = xp(find(f_xp==f_gbest_xp),:);
        if size(gbest_xp,1) > 1
            gbest_xp = gbest_xp(randperm(size(gbest_xp,1),1),:);
        end
        f_gbest_xp
        GBA_record3(iter)= f_gbest_xp;
        fprintf('Excecuting CCSO，the %dth iteration finished,%d pixels remained\n',iter,Pixels);
        if iter>200 && (GBA_record3(iter-200)-GBA_record3(iter))<0.001
            iter_finish = 1;
        end
    end
    for j = 1:size(xp,1)
        binarycode = transformRCToBinary(row,col,P,xp(j,:));
        realcode = find(binarycode==1);
        indicies = [indicies realcode];
    end

    %% Step3：Perform Fcls, remove the successfully unmixed pixels, and continue to unmixe the unsuccessfully unmixed pixels
    gbest_xp = transformRCToBinary(row,col,P,gbest_xp);
    endmember = img2d(:,gbest_xp==1);
    abundance = hyperFcls(img2d_reduced,endmember);
    remixed = endmember * abundance;
    residual_error = img2d_reduced - remixed;
    rmse = sqrt(sum(residual_error.^2)/L);  % RMSE of each pixel
    img2d_reduced(:,rmse<r) = [];
    [~,Pixels] = size(img2d_reduced);
    record{k+1,1} = GBA_record3;
    record{k+1,2} = rmse;
end







