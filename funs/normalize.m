function X_norm = normalize(X)
    % 全局归一化：将矩阵X的所有元素缩放到[0,1]
    % 输入：X为任意维度矩阵
    % 输出：X_norm为归一化后的矩阵
    
    min_val = min(X(:));  % 矩阵全局最小值
    max_val = max(X(:));  % 矩阵全局最大值
    
    if max_val == min_val  % 避免除以零（所有元素相同）
        X_norm = zeros(size(X));  % 或全设为0.5
    else
        X_norm = (X - min_val) / (max_val - min_val);
    end
end