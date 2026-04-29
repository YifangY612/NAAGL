function [UU,V,A,W,Z,E,iter,obj,y1] = naagl(X,Y,numanchor,alpha)
% m      : the number of anchor. the size of Z is m*n.
% X{p}   : n*di

%% initialize
maxIter = 30 ; % the number of iterations
IterMax = 50;

m = numanchor;
numclass = length(unique(Y));
numview = length(X);
numsample = size(Y,1);

Z = zeros(m,numsample); % m  * n
XX = [];
for p = 1 : numview
    X{p} = mapstd(X{p}',0,1);  % 0-1标准化  di*n
    XX = [XX;X{p}];
end
[XU,~,~]=svds(XX',m);  % 得到m个最大特征值对应的左奇异向量
rand('twister',12);
[IDX,~] = kmeans(XU,m, 'MaxIter',100,'Replicates',10);
for i = 1:numsample
    Z(IDX(i),i) = 1;
end

for i = 1:numview
   di = size(X{i},2); 
   E{i} = Z;
   A{i} = zeros(di,m);
   W{i} = eye(m);
end

flag = 1;
iter = 0;
%%
while flag
    iter = iter + 1;
    
    %% optimize Ai
    parfor iv=1:numview
        C = X{iv}*(W{iv}*Z+E{iv})';      
        [U,~,V] = svd(C,'econ');
        A{iv} = U*V';
    end


    %% optimize Wi
    parfor iv=1:numview
        R =Z*(X{iv}-A{iv}*E{iv})'*A{iv};      
        [U1,~,V1] = svd(R','econ');
        W{iv} = U1*V1';
    end
    
    %% optimize Ei
    parfor iv=1:numview
        H = (X{iv}-A{iv}*W{iv}*Z)'*A{iv};      
        for ii=1:numsample
            ut = H(ii,:)./(1+alpha);
            E{iv}(:,ii) = EProjSimplex_new(ut');
        end
    end
   
    %% optimize Z
    B = zeros(numsample,m);
    %G = zeros(m,numsample);
    for iv=1:numview
        B = B+(X{iv}-A{iv}*E{iv})'*A{iv}*W{iv};
        %G = G+E{iv}./numview;
    end
    B = B./numview;
    [P,~,~,y1] = coclustering_bipartite_fast_re(B, Z', numclass,IterMax);
    Z = P';
    
    term1 = 0;
    term2 = 0;
    for iv = 1:numview
        term1 = term1 + norm(X{iv}-A{iv}*(W{iv}*Z+E{iv}),'fro')^2;
        term2 = term2 + norm(E{iv},'fro')^2;
    end
    
    obj(iter) = term1+alpha*term2;
    
     
     % if iter>= maxIter
     %    [UU,~,V]=svd(Z','econ');
     %    flag = 0;
     % end
     
     if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-4 || iter>maxIter || obj(iter) < 1e-10)
        [UU,~,V]=svd(Z','econ');
        flag = 0;
     end



end
         
    
    % if iter<maxIter
    %        obj(iter+1:maxIter) =  obj(iter); 
    %    else    
    %        obj(maxIter+1:end) = [];
    % end
