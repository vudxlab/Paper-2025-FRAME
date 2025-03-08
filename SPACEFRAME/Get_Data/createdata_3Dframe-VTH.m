% Lexuanthang.official@gmail.com,Lexuanthang.official@outlook.com
% Le Xuan Thang, 2023
% Tutorial: dynamic analysis: direct method: frequency domain
% Units: m, N

%% Bước 1: Nhập mô hình
spaceframe;

% Nhập trường hợp hư hỏng 1 phần tử/ 1 lần
% 
% Matrix_Case = [0:9; % Thứ tự trường hợp / label
%     1 1 2 3 4 5 301 302 303 304  ; % phần tử chịu hư hỏng
%     0 10 20 30 40 50 40 30 20 10 ]; % Phần trăm hư hỏng

% Matrix_Case = [0:3; % Thứ tự trường hợp / label
%     1 1 2 3   ; % phần tử chịu hư hỏng
%     0 10 20 30 ]; % Phần trăm hư hỏng

Matrix_Case = [0:26; % Thứ tự trường hợp / label
    1 1 2 3 4 5 6 101 102 103 104 105 106 201 202 203 204 205 206 301 302 303 304 305 306 307 308; % phần tử chịu hư hỏng
    0 10 25 20 25 30 35 10 25 20 25 30 35 10 25 20 25 30 35 10 25 20 25 30 35 40 45]; % Phần trăm hư hỏng

Materials0 = Materials;
Elements0 = Elements;

for i = 1:size(Matrix_Case,2)
    Case = Matrix_Case(1, i);
    Element = Matrix_Case(2, i);
    Damage = Matrix_Case(3, i);
    Materials = Materials0;
    Elements = Elements0;
    % Update Materials and Elements based on the current case
    Materials = [Materials; 
                 3 35e9*(1-Damage/100) 0.2 2500 1.5e11];
    Elements(Elements(:,1) == Element, 4) = 3;

    % check frequency
    % Assembly of stiffness matrix K
    [K,M]=asmkm(Nodes,Elements,Types,Sections,Materials,DOF);

    % Eigenvalue problem
    nMode=12;
    [~,omega]=eigfem(K,M,nMode);
    frequency = omega/2/pi;

%% Gán tải
    %% BƯớc 2: Chọn loại phương tiện tác động lên cầu, chọn làn xe chạy, chọn tần số lấy mẫu
% Train/Vehicles phương tiện
%     P = [ -35 -145 -145;      % P1 P2 P3 / Lực trục
%     0 4.3 4.3];               % 0 l2 l3 khoảng cách giữa các trục
    seldof = [106.02,106.01];            % vị trí tác động 
    Pload =[100,100];
    P=nodalvalues(DOF,seldof,Pload);
    
%     DTBB =30;                 % [m] Distance train/truck to bridge before (front axle) khoảng cách giữa bánh xe và đầu cầu
%     V = 80*1000/3600;         % km/h --> [m/s] Velocity / Vận tốc chạy
%     LT = sum(P(2,:));         % [m] Length of train/truck 
% 
%     seldof = reprow([1500],1,8,[2])+0.03; %% *** Chọn làn mà xe chạy ***
% 
%     seldof = seldof(:);
%    dt = 0.002;            % Time step/resolution *** Chọn bước thời gian / tần số lấy mẫu***
    T=100;                   % chu kỳ lấy mẫu 10 phút
    dt=0.01;               % Tần số lấy mẫu
    N=T/dt;                 % số mẫu
    t = (0:N)*dt;           % Time axis (samples)
    % Sampling parameters: frequency domain
    F=1/dt;                 % Sampling frequency [Hz]
    df=1/T;                 % Frequency resolution
    f=(0:fix(N/2)-1)*df;    % Positive frequencies corresponding to FFT [Hz]
    Omega=2*pi*f;           % [rad/s] excitation frequency content

    % Excitation: transfer PLoad vector to nodalforce vector
    
    % (seldof --> all DOF)
%% Bước 3: Chọn các thông số lặp cho hàm gán lực vào kết cấu
%     startInterval = 2;        % Thời gian tàu bắt đầu vào cầu/ thời gian bắt đầu ghi dữ liệu 
%     nloop = 8;                % Số lần chạy của phương tiện trên cầu
%     Pulse = -500;             % Lực xung kích Pulse tác động lên cầu / Kết cấu
%     gap = 10;                 % Khoảng cách giữa các lần chạy của phương tiện trên cầu
%     PLoad = trainload(P,L,DTBB,V,dt,seldof,Nodes,startInterval,Pulse,T, nloop, gap,f0); % [seldof x N samples]

%%
    % Eigenvalue analysis
    nMode = 10;                     % Number of modes to take into account
    [phi,omega]=eigfem(K,M,nMode);  % Calculate eigenmodes and eigenfrequencies
    xi=0.02;                       % Constant modal damping ratio

    bm=phi.'*P;
    q=zeros(1,N);                   % Time history (1 * N)
%% sửa số lần gõ búa ở đây
%     q((t>=1.50) & (t<1.60) | (t>=10 & t <10.1) )=1;      % Time history (1 * N)
    q(mod(t, 10) == 0 & (t>0))=1;
%%---------------
    Q=fft(q);                       % Frequency content (1 * N)
    Q=Q(1:N/2);                     % Frequency content, positive freq (1 * N/2)
    Pm=bm*Q;

%     Pnodal = zeros(size(DOF,1),N);
%     for itime = 1:N
%         Pnodal(:,itime) = nodalvalues(DOF,seldof,PLoad(:,itime));
%     end
% 
%     % Modal excitation
%     Pm_ = phi.'*Pnodal; % [DOF,nMode] x [DOF, N samples]
% 
%     % Transfer nodal force vector time history to frequency domain
%     Q = zeros(nMode,fix(N/2)); % keep positive frequency ONLY
%     for inMode = 1:nMode
%         temp = fft(Pm_(inMode,:));
%         Q(inMode,:) = temp(1:fix(N/2));
%     end
% 
%     Pm = Q;

    % Modal analysis: calc. the modal transfer functions and the modal disp.
    [X,H]=msupf(omega,xi,Omega,Pm);     % Modal response, positive freq (nMode * N/2)
    % F-dom -> t-dom [inverse Fourier transform]
    X = [X, zeros(nMode,1), conj(X(:,end:-1:2))];
    x = ifft(X,[],2) ;                  % Modal response (nMode * N)
    x = x + abs(randn(size(x,1),size(x,2))).*(max(max(x))/1000);
    % Modal displacements -> nodal displacements
    u= phi*x ;                          % Nodal response (nDOF * N)

%% Bước 4: Chọn sensors/ điểm nodes xuất kết quả
%     
%    d = [reprow([101,201],1,1,[1,1])+0.03;reprow([101,201],1,1,[1,1])+0.02 ];
    d = [reprow([3,4,5,6],1,2,[100,100,100,100])+0.01; reprow([3,4,5,6],1,2,[100,100,100,100])+0.02];
%     d= [103.02]
    c = selectdof(DOF,d(:));
    u_c = c*u;

    acceleration = zeros(size(u_c));
%     for i_acc = 1:size(u_c/2,1)
    for i_acc = 1:size(u_c/2,1)
        acceleration(i_acc,:) = displacementToAcceleration(u_c(i_acc,:), dt);
%         figure;
%         plot(acceleration(i_acc,:));
%         title(["acc" i_acc]);
%         xlabel("samples");
%         ylabel("acc");
    end
%% Xem lực tác động

% % Movie
% figure;
% animdisp(Nodes,Elements,Types,DOF,u);
% animdisp(Nodes,Elements,Types,DOF,u(:,2));
    cur_dir = pwd();
    cd('D:\OneDrive\DXLaboratory\Papers\Trong nuoc\2024\Paper-2024-1DCNNLSTM-Khung\Data')
    filename = sprintf('spaceframe%d.mat', Case);
    save(filename, 'acceleration');
    cd(cur_dir)

% for i = 1:2
%     figure;
%     animdisp(Nodes,Elements,Types,DOF,u,'DispMax','off','LineWidth',2)
%     title({['Hình thái dao động ' num2str(i)], ['f= ' num2str(f0(i)) ' Hz']})
% end
end
