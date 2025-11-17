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

% --- 碱液 (T_f) 物性参数 ---
lambda_f = 0.68; % 碱液热导率 (W/m·K)
rho_f = 1200; % 碱液密度 (kg/m³)
cp_f = 3800; % 碱液比热容 (J/kg·K)
channel_depth_f = 2e-3; % 碱液流道深度 (m)

% --- !!! 新增: 碱液动力粘度 (Pa·s) (!! 这是一个假设值, 请修改 !!) !! ---
mu_f = 1e-3; % 例如: 1e-3 (Pa·s)

% --- 碱液流量 ---
Q_f_Lh = 500; % 碱液流量 (L/h)
Q_f = Q_f_Lh * (0.001 / 3600); % 转换为 m³/s (0.001 m³/L, 3600 s/h)

% --- !!! 新增: 计算 Pr (常数) !!! ---
Pr_f = (mu_f * cp_f) / lambda_f;
fprintf('计算得到的普朗特数 Pr = %f\n', Pr_f);

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 2. 派生参数和几何函数定义
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
R = D / 2;
H = D; % 总高度 (z从0到H)

% --- 几何匿名函数 (z 为高度) ---
L = @(z) 2 * sqrt(R^2 - (R-z).^2);
Ac_bpp = @(z) L(z) * T_thickness;
P_f_func = @(z) 2 * L(z);
P_air_func = @(z) 2 * T_thickness * ones(size(z));
A_f_func = @(z) L(z) * channel_depth_f;

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
P_f_nodes = P_f_func(z_nodes);
P_air_nodes = P_air_func(z_nodes);

% --- 在节点界面计算 ---
Ac_bpp_faces = Ac_bpp(z_faces);
A_f_faces = A_f_func(z_faces);
D_T_faces = lambda .* Ac_bpp_faces ./ dz;
D_f_faces = lambda_f .* A_f_faces ./ dz;
F_dot_f = rho_f * cp_f * Q_f;
F_f_faces = F_dot_f * ones(size(z_faces));

% --- !!! 新增: 计算 h_f(z) (在节点中心) !!! ---
% 特征长度 l = z
l_nodes = z_nodes;
% l_nodes = 1.9/100;

% 1. 碱液流通截面积 (在节点)
A_f_nodes = A_f_func(z_nodes);
% 2. 碱液流速 (在节点)
v_f_nodes = Q_f ./ A_f_nodes;
% 3. 雷诺数 Re_l (l=z)
Re_f_nodes = (rho_f * v_f_nodes .* l_nodes) / mu_f;
% 4. 努塞尔特数 Nu_l
Nu_f_nodes = 0.664 * (Re_f_nodes .^ 0.5) .* (Pr_f ^ (1/3));
% 5. 对流系数 h_f(z)
h_f_nodes = (Nu_f_nodes * lambda_f) ./ l_nodes;

% 6. 处理 z=0 (第一个节点) 的潜在问题
if any(isnan(h_f_nodes)) || any(isinf(h_f_nodes))
    fprintf('警告: h_f 计算在 z=0 附近出现 Inf/NaN。\n');
    % 简单的修正：将 NaN/Inf 替换为下一个有效值
    h_f_nodes(isnan(h_f_nodes)) = h_f_nodes(find(~isnan(h_f_nodes), 1, 'first'));
    h_f_nodes(isinf(h_f_nodes)) = h_f_nodes(find(~isinf(h_f_nodes), 1, 'first'));
end

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 5. 初始化
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
T = ones(N, 1) * t_air; % BPP 温度场
T_f = ones(N, 1) * T_bottom_fixed; % 碱液温度场
aW_T = zeros(N, 1); aE_T = zeros(N, 1); aP_T = zeros(N, 1);
SP_T = zeros(N, 1); SU_T = zeros(N, 1);
aW_f = zeros(N, 1); aE_f = zeros(N, 1); aP_f = zeros(N, 1);
SP_f = zeros(N, 1); SU_f = zeros(N, 1);
fprintf('开始 FVM 耦合迭代...\n');

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 6. 外循环 (耦合求解)
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
for iter_outer = 1:max_outer_iter

    T_old = T;
    T_f_old = T_f;
    
    % --- 6A. 组装并求解 BPP 温度 T ---
    
    % (a) 内部节点系数
    for n = 2:N-1
        aW_T(n) = D_T_faces(n); 
        aE_T(n) = D_T_faces(n+1);
    end
    
    % (b) 计算源项 (SU, SP) - 依赖 T_f_old
    % S_T = (q_dot*Ac_bpp) + h_f*P_f*T_f + h_air*P_air*t_air
    %       - (h_f*P_f + h_air*P_air) * T
    SU_T = (q_dot_vol .* Ac_bpp_nodes + ...
            h_f_nodes .* P_f_nodes .* T_f_old + ...  % <-- 已修改
            h_air .* P_air_nodes .* t_air) * dz;
    SP_T = -(h_f_nodes .* P_f_nodes + h_air .* P_air_nodes) * dz; % <-- 已修改
    
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

    % --- 6B. 组装并求解 碱液 温度 T_f ---

    % (a) 内部节点系数 (UDS)
    for n = 2:N-1
        aW_f(n) = D_f_faces(n) + F_f_faces(n);
        aE_f(n) = D_f_faces(n+1);
    end

    % (b) 计算源项 (SU, SP) - 依赖 *新计算的 T*
    % S_f = (h_f * P_f * T) - (h_f * P_f) * T_f
    SU_f = (h_f_nodes .* P_f_nodes .* T) * dz; % <-- 已修改
    SP_f = -(h_f_nodes .* P_f_nodes) * dz; % <-- 已修改

    % (c) 计算 aP
    aP_f = aW_f + aE_f - SP_f;

    % (d) T_f 的边界条件
    % 节点 1 (Dirichlet: T_f = T_bottom_fixed)
    aW_f(1) = 0;
    aE_f(1) = 0;
    aP_f(1) = 1;
    SU_f(1) = T_bottom_fixed;
    SP_f(1) = 0;

    % 节点 N (Neumann: dTf/dz=0)
    aW_f(N) = D_f_faces(N) + F_f_faces(N);
    aE_f(N) = 0;
    aP_f(N) = aW_f(N) - SP_f(N);

    % (e) 内循环: 迭代求解 T_f
    for iter_inner = 1:max_inner_iter
        T_f_inner_old = T_f;
        T_f(1) = SU_f(1) / aP_f(1); % 固定为 80
        for n = 2:N-1
            T_f(n) = (aW_f(n) * T_f(n-1) + aE_f(n) * T_f(n+1) + SU_f(n)) / aP_f(n);
        end
        T_f(N) = (aW_f(N) * T_f(N-1) + SU_f(N)) / aP_f(N);
        if max(abs(T_f - T_f_inner_old)) < inner_tol
            break;
        end
    end
    
    % --- 6C. 检查外循环收敛 ---
    err_T = max(abs(T - T_old));
    err_Tf = max(abs(T_f - T_f_old));
    residual = max(err_T, err_Tf);
    
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
plot(z_nodes, T_f, 'r--', 'LineWidth', 2);
xlabel('高度 z (m)');
ylabel('温度 (°C)');
title('FVM 求解 BPP 和碱液温度分布');
grid on;

% % --- 子图 2: 计算得到的 h_f 分布 ---
% subplot(2, 1, 2); % 选择第2个
% plot(z_nodes, h_f_nodes, 'm-', 'LineWidth', 2);
% xlabel('高度 z (m)');
% ylabel('对流系数 h_{f} (W/m²·K)');
% title('计算得到的 h_{f} 沿高度的分布');
% grid on;
