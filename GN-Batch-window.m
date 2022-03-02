% AER1513: Assignment 3 - Q.5
clear;

% loading the given data and functions
load dataset3;
addpath('barfoot_tro14');

% defining limits
k1 = 1215;
k2 = 1714;

% counting the time-steps with visible landmarks
count = 0;
for k=k1:k2
    for j=1:size(y_k_j,3)
        if y_k_j(1,k,j) ~= -1
            count = count + 1;
        end
    end
end

% saving required sizes to define matrices
big_size = count * 4;
timesteps = k2-k1+1;
last = 6*timesteps;

% initialising rotation, position and transformation matrices
C_vk_i = zeros(3,3,timesteps);
r_vk_i = zeros(3,timesteps);
T_vk_i = zeros(4,4,timesteps);
T_vk_i(4,4,:) = 1;

% defining inputs and variances
varpi_vk_vk_i = [-1*v_vk_vk_i; -1*w_vk_vk_i];
vw_var = [v_var; w_var];

% defining transformtion from vehicle frame to camera frame
T_c_v = eye(4);
T_c_v(1:3,1:3) = C_c_v;
T_c_v(1:3,4) = -1*C_c_v*rho_v_c_v;

% initialising matrices for Batch estimation
F = zeros(6,6,timesteps);
H = zeros(big_size+last,last);
e_k = zeros(big_size+last,1);
R_k = diag(y_var);
W = zeros(big_size+last);
i_R = last + 1;
del_x = ones(last,1);
Q_k = zeros(last);

% camera model
M = [fu 0 cu 0; 0 fv cv 0; fu 0 cu -fu*b; 0 fv cv 0];

% initiation
for i = k1:k2
    
    % current index for defined matrices
    a = i-k1+1;
    
    % motion model calculations
    if a == 1
        % at operating point
        C_vk_i(:,:,a) = vec2rot(-1*theta_vk_i(:,i));
        r_vk_i(:,a) = r_i_vk_i(:,i);
        T_vk_i(1:3,1:3,a) = C_vk_i(:,:,a);
        T_vk_i(1:3,4,a) = -1*C_vk_i(:,:,a)*r_vk_i(:,a);
    else
        d_vk = v_vk_vk_i(:,i-1)*(t(i)-t(i-1));
        psi = w_vk_vk_i(:,i-1)*(t(i)-t(i-1));
        angle = norm(psi);
        axis = psi/angle;
        psi_vk = cos(angle)*eye(3) + (1-cos(angle))*(axis*axis') - sin(angle)*hat(axis);
        r_vk_i(:,a) = r_vk_i(:,a-1) + (C_vk_i(:,:,a-1)')*d_vk;
        C_vk_i(:,:,a) = psi_vk*C_vk_i(:,:,a-1);
        T_vk_i(1:3,1:3,a) = C_vk_i(:,:,a);
        T_vk_i(1:3,4,a) = -1*C_vk_i(:,:,a)*r_vk_i(:,a);
        Xi_k = vec2tran((t(i)-t(i-1))*varpi_vk_vk_i(:,i-1));
        e_k(6*a-5:6*a) = round(tran2vec(Xi_k*T_vk_i(:,:,a-1)*tranInv(T_vk_i(:,:,a))));
        F(:,:,a-1) = tranAd(T_vk_i(:,:,a)*tranInv(T_vk_i(:,:,a-1)));
        H(6*a-5:6*a, 6*a-11:6*a-6) = -1*F(:,:,a-1);
    end
    
    H(6*a-5:6*a,6*a-5:6*a) = eye(6);
    Q(6*a-5:6*a,6*a-5:6*a) = ((t(i)-t(i-1))^2)*diag(vw_var);
    W(6*a-5:6*a,6*a-5:6*a) = Q(6*a-5:6*a,6*a-5:6*a);
    flag = 0;
    
    % measurement model calculations
    for j = 1:20
        % considering only visible landmarks
        if y_k_j(1,i,j) ~= -1
%             p_ck_pj_ck = C_c_v*(vec2rot(-theta_vk_i(:,i))*(rho_i_pj_i(:,j) - r_i_vk_i(:,i))-rho_v_c_v);
            rho = T_c_v * T_vk_i(:,:,a) * [rho_i_pj_i(:,j);1];
            S_j = camJac(M,rho);
            Z_j_ = (T_vk_i(:,:,a)*[rho_i_pj_i(:,j); 1]);
            Z_j = [Z_j_(4)*eye(3), -1*hat(Z_j_(1:3)); zeros(1,6)];
            G_j = S_j*T_c_v*Z_j;
            y_k_j_ = M*(1/rho(3))*rho;
            e_y = y_k_j(:,i,j) - y_k_j_;
            if flag == 0
                G_j_k = G_j;
                e_y_k = e_y;
                flag = 1;
            else
                G_j_k = [G_j_k; G_j];
                e_y_k = [e_y_k; e_y];
            end
            W(i_R:i_R+3,i_R:i_R+3) = R_k;
            i_R = i_R + 4;
        end
    end
    
    if flag == 1
        H(last+1:last+size(G_j_k,1),6*a-5:6*a) = G_j_k;    
        e_k(last+1:last+size(e_y_k,1)) = e_y_k;
        last = last + size(G_j_k,1);
    end
    
end

% Gauss-Newton
count = 0;
check = 1;
while check > 0.01

    count = count + 1;
    last = 6*timesteps;

    A = H' / W * H;
    B = H' / W * real(e_k);
    del_x = A \ B;
    check = norm(del_x);
    
    for l = 1:6:size(del_x,1)
        a = (l+5)/6;
        i = a + k1 - 1;
        T_vk_i(:,:,a) = vec2tran(del_x(l:l+5))*T_vk_i(:,:,a);
    
        if a > 1
            e_k(6*a-5:6*a) = round(tran2vec(Xi_k*T_vk_i(:,:,a-1)*tranInv(T_vk_i(:,:,a))));
            F(:,:,a-1) = tranAd(T_vk_i(:,:,a)*tranInv(T_vk_i(:,:,a-1)));
            H(6*a-5:6*a, 6*a-11:6*a-6) = -1*F(:,:,a-1);
        end
        
        % measurement model calculations
        flag = 0;
        for j = 1:20        
            % considering only visible landmarks
            if y_k_j(1,i,j) ~= -1
%                 p_ck_pj_ck = C_c_v*(vec2rot(-theta_vk_i(:,i))*(rho_i_pj_i(:,j) - r_i_vk_i(:,i))-rho_v_c_v);
                rho = T_c_v * T_vk_i(:,:,a) * [rho_i_pj_i(:,j);1];
                S_j = camJac(M,rho);
                Z_j_ = (T_vk_i(:,:,a)*[rho_i_pj_i(:,j); 1]);
                Z_j = [Z_j_(4)*eye(3), -1*hat(Z_j_(1:3)); zeros(1,6)];
                G_j = S_j*T_c_v*Z_j;
                y_k_j_ = M*(1/rho(3))*rho;
                e_y = y_k_j(:,i,j) - y_k_j_;
                if flag == 0
                    G_j_k = G_j;
                    e_y_k = e_y;
                    flag = 1;
                else
                    G_j_k = [G_j_k; G_j];
                    e_y_k = [e_y_k; e_y];
                end
            end
        end
        
        if flag == 1
            H(last+1:last+size(G_j_k,1),6*a-5:6*a) = G_j_k;    
            e_k(last+1:last+size(e_y_k,1)) = e_y_k;
            last = last + size(G_j_k,1);
        end
    end
end

% extracting rotation and position matrices
dC_k = zeros(3,3,timesteps);
del_r_k = zeros(3,timesteps);
del_theta_k = zeros(3,timesteps);
del_ = zeros(6,timesteps);
C_vk_i = T_vk_i(1:3,1:3,:);
for k=1:timesteps
    r_vk_i = -1 * T_vk_i(1:3,1:3,k)' * T_vk_i(1:3,4,k);
    del_r_k(:,k) = r_vk_i - r_i_vk_i(:,k+k1-1);
    dC_k(:,:,k) = C_vk_i(:,:,k) * vec2rot(theta_vk_i(:,k+k1-1))';
    diff = eye(3) - dC_k(:,:,k);
    del_theta_k(:,k) = [diff(3,2); diff(1,3); diff(2,1)];
end
del_k = [del_r_k; del_theta_k];            

% % plotting dead-reckoned values from k1 to k2, with ground truth
% figure(1)
% hold on;
% plot3(rho_i_pj_i(1,:), rho_i_pj_i(2,:), rho_i_pj_i(3,:), 'r.', 'markersize', 15);
% plot3(r_vk_i(1,:), r_vk_i(2,:), r_vk_i(3,:), 'b');
% plot3(r_i_vk_i(1,k1:k2), r_i_vk_i(2,k1:k2), r_i_vk_i(3,k1:k2), 'g');
% grid on;

% standard deviations
P = inv(A);
std_k = sqrt(diag(P));
std_k_1 = zeros(1,timesteps);
std_k_2 = zeros(1,timesteps);
std_k_3 = zeros(1,timesteps);
std_k_4 = zeros(1,timesteps);
std_k_5 = zeros(1,timesteps);
std_k_6 = zeros(1,timesteps);

for p = 1:timesteps
    std_k_1(p) = std_k(6*p-5);
    std_k_2(p) = std_k(6*p-4);
    std_k_3(p) = std_k(6*p-3);
    std_k_4(p) = std_k(6*p-2);
    std_k_5(p) = std_k(6*p-1);
    std_k_6(p) = std_k(6*p);
end

% plotting it all
figure(2)
plot(del_r_k(1,:), 'linewidth', 2);
title("error in linear position in x-direction v/s time", 'fontsize', 20);
ylabel("\deltar_x", 'fontsize', 20);
xlabel("time-steps", 'fontsize', 20);
hold on;
plot(3*std_k_1, ':', 'linewidth', 2);
plot(-3*std_k_1, ':', 'linewidth', 2);
legend("\deltar_x", "+3\sigma_r_x", "-3\sigma_r_x", 'fontsize', 20);
grid on;

figure(3)
plot(del_r_k(2,:), 'linewidth', 2);
title("error in linear position in y-direction v/s time", 'fontsize', 20);
ylabel("\deltar_y", 'fontsize', 20);
xlabel("time-steps", 'fontsize', 20);
hold on;
plot(3*std_k_2, ':', 'linewidth', 2);
plot(-3*std_k_2, ':', 'linewidth', 2);
legend("\deltar_y", "+3\sigma_r_y", "-3\sigma_r_y", 'fontsize', 20);
grid on;

figure(4)
plot(del_r_k(3,:), 'linewidth', 2);
title("error in linear position in z-direction v/s time", 'fontsize', 20);
ylabel("\deltar_z", 'fontsize', 20);
xlabel("time-steps", 'fontsize', 20);
hold on;
plot(3*std_k_3, ':', 'linewidth', 2);
plot(-3*std_k_3, ':', 'linewidth', 2);
legend("\deltar_z", "+3\sigma_r_z", "-3\sigma_r_z", 'fontsize', 20);
grid on;

figure(5)
plot(del_theta_k(1,:), 'linewidth', 2);
title("error in angular position about x-direction v/s time", 'fontsize', 20);
ylabel("\delta\theta_x", 'fontsize', 20);
xlabel("time-steps", 'fontsize', 20);
hold on;
plot(3*std_k_4, ':', 'linewidth', 2);
plot(-3*std_k_4, ':', 'linewidth', 2);
legend("\delta\theta_x", "+3\sigma_\theta_x", "-3\sigma_\theta_x", 'fontsize', 20);
grid on;

figure(6)
plot(del_theta_k(2,:), 'linewidth', 2);
title("error in angular position about y-direction v/s time", 'fontsize', 20);
ylabel("\delta\theta_y", 'fontsize', 20);
xlabel("time-steps", 'fontsize', 20);
hold on;
plot(3*std_k_5, ':', 'linewidth', 2);
plot(-3*std_k_5, ':', 'linewidth', 2);
legend("\delta\theta_y", "+3\sigma_\theta_y", "-3\sigma_\theta_y", 'fontsize', 20);
grid on;

figure(7)
plot(del_theta_k(3,:), 'linewidth', 2);
title("error in angular position about z-direction v/s time", 'fontsize', 20);
ylabel("\delta\theta_z", 'fontsize', 20);
xlabel("time-steps", 'fontsize', 20);
hold on;
plot(3*std_k_6, ':', 'linewidth', 2);
plot(-3*std_k_6, ':', 'linewidth', 2);
legend("\delta\theta_z", "+3\sigma_\theta_z", "-3\sigma_\theta_z", 'fontsize', 20);
grid on;

function [T_inv] = tranInv(T)
    C = T(1:3, 1:3);
    r = T(1:3, 4);
    T_inv = [C' -C'*r; 0 0 0 1];
end