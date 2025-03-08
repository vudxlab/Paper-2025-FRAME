% StaBIL manual
% Example 2.1: dynamic analysis: direct method: frequency domain
% Units: m, N

% Assembly of M and K
tutorialdyna;

% Sampling parameters
N=2048;                         % Number of samples
dt=0.002;                       % Time step
T=N*dt;                         % Period
F=N/T;                          % Sampling frequency
df=1/T;                         % Frequency resolution
t=[0:N-1]*dt;                   % Time axis
f=[0:N/2-1]*df;                 % Positive frequencies corresponding to FFT [Hz]
Omega=2*pi*f;                   % Idem [rad/s]

% Eigenvalue analysis
nMode=12;                       % Number of modes to take into account
[phi,omega]=eigfem(K,M,nMode);  % Calculate eigenmodes and eigenfrequencies
xi=0.07;                        % Constant modal damping ratio

% Excitation
bm=phi.'*b;                     % Spatial distribution, modal (nMode * 1)
q=zeros(1,N);                   % Time history (1 * N)
q((t>=0.50) & (t<0.60))=1;      % Time history (1 * N)
Q=fft(q);                       % Frequency content (1 * N)
Q=Q(1:N/2);                     % Frequency content, positive freq (1 * N/2)
Pm=bm*Q;                        % Modal excitation, positive freq (nMode * N/2)

% Modal analysis
[X,H]=msupf(omega,xi,Omega,Pm); % Modal response, positive freq (nMode * N/2)

% F-dom -> t-dom
X=[X, zeros(nMode,1), conj(X(:,end:-1:2))];
x=ifft(X,[],2);                 % Modal response (nMode * N)

% Modal displacements -> nodal displacements
u=phi*x;                        % Nodal response (nDOF * N)

% Figures
figure;
subplot(3,2,1);
plot(t,q,'.-');
xlim([0 4.1])
ylim([0 1.2]);
title('Excitation time history');
xlabel('Time [s]');
ylabel('Force [N/m]');

subplot(3,2,2);
plot(f,abs(Q)/F,'.-');
title('Excitation frequency content');
xlabel('Frequency [Hz]');
ylabel('Force [N/m/Hz]');

subplot(3,2,4);
plot(f,abs(H),'.-');
title('Modal transfer function');
xlabel('Frequency [Hz]');
ylabel('Displacement [m/N]');
legend([repmat('Mode ',nMode,1) num2str([1:nMode].')]);

subplot(3,2,6);
plot(f,abs(X(:,1:N/2))/F,'.-');
title('Modal response');
xlabel('Frequency [Hz]');
ylabel('Displacement [m kg^{0.5}/Hz]');

subplot(3,2,5);
plot(t,x);
title('Modal response (calculation in f-dom)');
xlabel('Time [s]');
xlim([0 4.1])
ylabel('Displacement [m kg^{0.5}]');

figure;
plot(t,x);
title('Modal response (calculation in f-dom)');
xlabel('Time [s]');
xlim([0 4.1])
ylabel('Displacement [m kg^{0.5}]');
legend([repmat('Mode ',nMode,1) num2str([1:nMode].')]);

figure;
c=selectdof(DOF,[9.01; 13.02; 17.02]);
plot(t,c*u);
title('Nodal response (calculation in f-dom)');
xlabel('Time [s]');
xlim([0 4.1])
ylabel('Displacement [m]');
legend('9.01','13.02','17.02');

% Movie
figure;
animdisp(Nodes,Elements,Types,DOF,u);

% Display
disp('Maximum modal response');
disp(max(abs(x),[],2));

disp('Maximum nodal response 9.01 13.02 17.02');
disp(max(abs(c*u),[],2));
