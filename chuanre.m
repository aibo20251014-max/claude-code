clear; clc; close all;

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 1. 物理和几何参数
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% --- BPP (T) 和环境参数 ---
t_air = 25; % 环境空气温度 (°C)
P_total = 7100; % BPP 总产热功率 (W)
lambda = 25; % BPP 热导率 (W/m·K)
% h_ele = 200; % (已删除) BPP 与碱液的对流系数 (W/m²·K)
h_air = 20; % BPP 与空气的对流系数 (W/m²·K)
D = 1.9; % BPP 高度/直径 (m)
T_thickness = 5e-3; % BPP 厚度 (m)
T_bottom_fixed = 80; % 底部 碱液 入口温度 (°C)

% --- 碱液 (T_ele) 物性参数 ---
lambda_ele = 0.68; % 碱液热导率 (W/m·K)
rho_ele = 1200; % 碱液密度 (kg/m³)
cp_ele = 3800; % 碱液比热容 (J/kg·K)
channel_depth_ele = 2e-3; % 碱液流道深度 (m)

% --- !!! 新增: 碱液动力粘度 (Pa·s) (!! 这是一个假设值, 请修改 !!) !! ---
mu_ele = 1e-3; % 例如: 1e-3 (Pa·s) 

% --- 碱液流量 ---
Q_ele_Lh = 500; % 碱液流量 (L/h)
Q_ele = Q_ele_Lh * (0.001 / 3600); % 转换为 m³/s (0.001 m³/L, 3600 s/h)

% --- !!! 新增: 计算 Pr (常数) !!! ---
Pr_ele = (mu_ele * cp_ele) / lambda_ele; 
fprintf('计算得到的普朗特数 Pr = %f\n', Pr_ele);

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 2. 派生参数和几何函数定义
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
R = D / 2;
H = D; % 总高度 (z从0到H)

% --- 几何匿名函数 (z 为高度) ---
L = @(z) 2 * sqrt(R^2 - (R-z).^2);
Ac_bpp = @(z) L(z) * T_thickness;
P_ele_func = @(z) 2 * L(z);
P_air_func = @(z) 2 * T_thickness * ones(size(z));
A_ele_func = @(z) L(z) * channel_depth_ele;

% --- 体积热源 ---
V_total_BPP = integral(Ac_bpp, 0, H);
q_dot_vol = P_total / V_total_BPP; % 体积热源 (W/m³)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 3. FVM 网格和迭代参数
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
N = 100; % 网格节点数
dz = H / N; % 网格间距
z_nodes = (dz/2 : dz : H-dz/2)';
z_faces = (0 : dz : H)';

% --- 迭代控制 ---
max_outer_iter = 500;
outer_tol = 1e-6;
max_inner_iter = 100;
inner_tol = 1e-8;

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 4. 预计算几何和物理参数
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% --- 在节点中心计算 ---
Ac_bpp_nodes = Ac_bpp(z_nodes);
P_ele_nodes = P_ele_func(z_nodes);
P_air_nodes = P_air_func(z_nodes);

% --- 在节点界面计算 ---
Ac_bpp_faces = Ac_bpp(z_faces);
A_ele_faces = A_ele_func(z_faces);
D_T_faces = lambda .* Ac_bpp_faces ./ dz;
D_ele_faces = lambda_ele .* A_ele_faces ./ dz;
F_dot_ele = rho_ele * cp_ele * Q_ele;
F_ele_faces = F_dot_ele * ones(size(z_faces));

% --- !!! 新增: 计算 h_ele(z) (在节点中心) !!! ---
% 特征长度 l = z
l_nodes = z_nodes;
% l_nodes = 1.9/100;

% 1. 碱液流通截面积 (在节点)
A_ele_nodes = A_ele_func(z_nodes);
% 2. 碱液流速 (在节点)
v_ele_nodes = Q_ele ./ A_ele_nodes;
% 3. 雷诺数 Re_l (l=z)
Re_ele_nodes = (rho_ele * v_ele_nodes .* l_nodes) / mu_ele;
% 4. 努塞尔特数 Nu_l
Nu_ele_nodes = 0.664 * (Re_ele_nodes .^ 0.5) .* (Pr_ele ^ (1/3));
% 5. 对流系数 h_ele(z)
h_ele_nodes = (Nu_ele_nodes * lambda_ele) ./ l_nodes;

% 6. 处理 z=0 (第一个节点) 的潜在问题
if any(isnan(h_ele_nodes)) || any(isinf(h_ele_nodes))
    fprintf('警告: h_ele 计算在 z=0 附近出现 Inf/NaN。\n');
    % 简单的修正：将 NaN/Inf 替换为下一个有效值
    h_ele_nodes(isnan(h_ele_nodes)) = h_ele_nodes(find(~isnan(h_ele_nodes), 1, 'first'));
    h_ele_nodes(isinf(h_ele_nodes)) = h_ele_nodes(find(~isinf(h_ele_nodes), 1, 'first'));
end

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 5. 初始化
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
T = ones(N, 1) * t_air; % BPP 温度场
T_ele = ones(N, 1) * T_bottom_fixed; % 碱液温度场
aW_T = zeros(N, 1); aE_T = zeros(N, 1); aP_T = zeros(N, 1);
SP_T = zeros(N, 1); SU_T = zeros(N, 1);
aW_ele = zeros(N, 1); aE_ele = zeros(N, 1); aP_ele = zeros(N, 1);
SP_ele = zeros(N, 1); SU_ele = zeros(N, 1);
fprintf('开始 FVM 耦合迭代...\n');

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 6. 外循环 (耦合求解)
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
for iter_outer = 1:max_outer_iter
    
    T_old = T;
    T_ele_old = T_ele;
    
    % --- 6A. 组装并求解 BPP 温度 T ---
    
    % (a) 内部节点系数
    for n = 2:N-1
        aW_T(n) = D_T_faces(n); 
        aE_T(n) = D_T_faces(n+1);
    end
    
    % (b) 计算源项 (SU, SP) - 依赖 T_ele_old
    % S_T = (q_dot*Ac_bpp) + h_ele*P_ele*T_ele + h_air*P_air*t_air
    %       - (h_ele*P_ele + h_air*P_air) * T
    SU_T = (q_dot_vol .* Ac_bpp_nodes + ...
            h_ele_nodes .* P_ele_nodes .* T_ele_old + ...  % <-- 已修改
            h_air .* P_air_nodes .* t_air) * dz;
    SP_T = -(h_ele_nodes .* P_ele_nodes + h_air .* P_air_nodes) * dz; % <-- 已修改
    
    % (c) 计算 aP
    aP_T = aW_T + aE_T - SP_T;
    
    % (d) T 的边界条件 (Neumann: dT/dz=0)
    aW_T(1) = 0; 
    aE_T(1) = D_T_faces(2);
    aP_T(1) = aE_T(1) - SP_T(1);
    aW_T(N) = D_T_faces(N);
    aE_T(N) = 0; 
    aP_T(N) = aW_T(N) - SP_T(N);

    % (e) 内循环: 迭代求解 T
    for iter_inner = 1:max_inner_iter
        T_inner_old = T;
        T(1) = (aE_T(1) * T(2) + SU_T(1)) / aP_T(1);
        for n = 2:N-1
            T(n) = (aW_T(n) * T(n-1) + aE_T(n) * T(n+1) + SU_T(n)) / aP_T(n);
        end
        T(N) = (aW_T(N) * T(N-1) + SU_T(N)) / aP_T(N);
        if max(abs(T - T_inner_old)) < inner_tol
            break;
        end
    end

    % --- 6B. 组装并求解 碱液 温度 T_ele ---
    
    % (a) 内部节点系数 (UDS)
    for n = 2:N-1
        aW_ele(n) = D_ele_faces(n) + F_ele_faces(n);
        aE_ele(n) = D_ele_faces(n+1);
    end
    
    % (b) 计算源项 (SU, SP) - 依赖 *新计算的 T*
    % S_ele = (h_ele * P_ele * T) - (h_ele * P_ele) * T_ele
    SU_ele = (h_ele_nodes .* P_ele_nodes .* T) * dz; % <-- 已修改
    SP_ele = -(h_ele_nodes .* P_ele_nodes) * dz; % <-- 已修改
    
    % (c) 计算 aP
    aP_ele = aW_ele + aE_ele - SP_ele;

    % (d) T_ele 的边界条件
    % 节点 1 (Dirichlet: T_ele = T_bottom_fixed)
    aW_ele(1) = 0;
    aE_ele(1) = 0;
    aP_ele(1) = 1;
    SU_ele(1) = T_bottom_fixed;
    SP_ele(1) = 0;
    
    % 节点 N (Neumann: dTele/dz=0)
    aW_ele(N) = D_ele_faces(N) + F_ele_faces(N);
    aE_ele(N) = 0;
    aP_ele(N) = aW_ele(N) - SP_ele(N);
    
    % (e) 内循环: 迭代求解 T_ele
    for iter_inner = 1:max_inner_iter
        T_ele_inner_old = T_ele;
        T_ele(1) = SU_ele(1) / aP_ele(1); % 固定为 80
        for n = 2:N-1
            T_ele(n) = (aW_ele(n) * T_ele(n-1) + aE_ele(n) * T_ele(n+1) + SU_ele(n)) / aP_ele(n);
        end
        T_ele(N) = (aW_ele(N) * T_ele(N-1) + SU_ele(N)) / aP_ele(N);
        if max(abs(T_ele - T_ele_inner_old)) < inner_tol
            break;
        end
    end
    
    % --- 6C. 检查外循环收敛 ---
    err_T = max(abs(T - T_old));
    err_Tele = max(abs(T_ele - T_ele_old));
    residual = max(err_T, err_Tele);
    
    if mod(iter_outer, 25) == 0
        fprintf('迭代 %d, 最大残差 = %g\n', iter_outer, residual);
    end
    
    if residual < outer_tol
        fprintf('耦合求解收敛于 %d 次迭代。\n', iter_outer);
        break;
    end
end
if iter_outer == max_outer_iter
    fprintf('警告: 达到最大迭代次数 %d，求解可能未收敛。\n', max_outer_iter);
end

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 7. 结果绘图
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
figure;

% --- 子图 1: 温度分布 ---
% subplot(2, 1, 1); % 分为2行1列，选择第1个
plot(z_nodes, T, 'b-', 'LineWidth', 2);
hold on;
plot(z_nodes, T_ele, 'r--', 'LineWidth', 2);
xlabel('高度 z (m)');
ylabel('温度 (°C)');
title('FVM 求解 BPP 和碱液温度分布');
grid on;

% % --- 子图 2: 计算得到的 h_ele 分布 ---
% subplot(2, 1, 2); % 选择第2个
% plot(z_nodes, h_ele_nodes, 'm-', 'LineWidth', 2);
% xlabel('高度 z (m)');
% ylabel('对流系数 h_{ele} (W/m²·K)');
% title('计算得到的 h_{ele} 沿高度的分布');
% grid on;
