
%% Bước 1: Nhập mô hình
spaceframe;

% Nhập trường hợp hư hỏng 1 phần tử/ 1 lần

Matrix_Case = [0:9; % Thứ tự trường hợp / label
    1 1 2 3 4 5 301 302 303 304  ; % phần tử chịu hư hỏng
    0 10 20 30 40 50 40 30 20 10 ]; % Phần trăm hư hỏng


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
    seldof = [106.02,106.01];            % vị trí tác động 
    Pload =[100,100];
    P=nodalvalues(DOF,seldof,Pload);
    

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
    % Eigenvalue analysis
    nMode = 10;                     % Number of modes to take into account
    [phi,omega]=eigfem(K,M,nMode);  % Calculate eigenmodes and eigenfrequencies
    xi=0.02;                       % Constant modal damping ratio

    bm=phi.'*P;
    q=zeros(1,N);                   % Time history (1 * N)
%% sửa số lần gõ búa ở đây
%     q((t>=1.50) & (t<1.60) | (t>=10 & t <10.1) )=1;      % Time history (1 * N)
    for a = 9:10
        q(mod(t, a) == 0 & (t>0))=1;
    %%---------------
        Q=fft(q);                       % Frequency content (1 * N)
        Q=Q(1:N/2);                     % Frequency content, positive freq (1 * N/2)
        Pm=bm*Q;
    
        % Modal analysis: calc. the modal transfer functions and the modal disp.
        [X,H]=msupf(omega,xi,Omega,Pm);     % Modal response, positive freq (nMode * N/2)
        % F-dom -> t-dom [inverse Fourier transform]
        X = [X, zeros(nMode,1), conj(X(:,end:-1:2))];
        x = ifft(X,[],2) ;                  % Modal response (nMode * N)
        x = x + abs(randn(size(x,1),size(x,2))).*(max(max(x))/500);
        % Modal displacements -> nodal displacements
        u= phi*x ;                          % Nodal response (nDOF * N)
    
    %% Bước 4: Chọn sensors/ điểm nodes xuất kết quả
    %     
    %    d = [reprow([101,201],1,1,[1,1])+0.03;reprow([101,201],1,1,[1,1])+0.02 ];
%         d = [reprow([3,4,5,6],1,2,[100,100,100,100])+0.01; reprow([3,4,5,6],1,2,[100,100,100,100])+0.02];
        d = [reprow([3,6],1,2,[100,100])+0.01; reprow([4,5],1,2,[100,100])+0.02];
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
%         cd('D:\OneDrive\DXLaboratory\Papers\Trong nuoc\2024\Paper-2024-1DCNNLSTM-Khung\Data')
        cd('C:\Users\vulv2\OneDrive\DXLaboratory\Papers\Trong nuoc\2024\Paper-2024-1DCNNLSTM-Khung\Code')
        filename = sprintf('spaceframe%d.mat', i);
        if exist(filename, 'file') == 2
            % Nếu tệp đã tồn tại, nối thêm vào
            existing_data = load(filename); % Tải dữ liệu từ tệp đã tồn tại
            data_cell = struct2cell(existing_data);
            existing_data = cell2mat(data_cell);
            acceleration = cat(1, existing_data, acceleration); % Nối dữ liệu mới vào dữ liệu đã có
            save(filename, 'acceleration'); % Lưu dữ liệu đã nối vào tệp
            disp('Dữ liệu đã được nối vào tệp .mat đã tồn tại.')
        else
            save(filename, 'acceleration');
            cd(cur_dir)
                
        end
        
    % for i = 1:2
    %     figure;
    %     animdisp(Nodes,Elements,Types,DOF,u,'DispMax','off','LineWidth',2)
    %     title({['Hình thái dao động ' num2str(i)], ['f= ' num2str(f0(i)) ' Hz']})
    % end
    end
end
