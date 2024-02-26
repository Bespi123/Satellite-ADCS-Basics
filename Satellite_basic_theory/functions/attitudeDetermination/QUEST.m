%function q = QUEST(fb, mb, fn, mn, wf, wm)
function q = QUEST(fb, mn, w)
% q = QUEST(Sensors, Parameters)
% Function implements QUEST algorithm using measurements
% from three-component accelerometer with orthogonal axes and vector
% magnetometer
%
%   Input arguments:
%   fb  - Vectors in body frame [3xn]
%   mn  - Vectors in navigation frame [3xn]
%   w   - Vector weight [1xn]

%   Output arguments:
%   Cbn - estimated Direction Cosines Matrix (DCM)


%% Wahba's problem 
%%
% We need to find $C_n^b$ matrix to minimize loss function:

%%
% $L(A) = \frac{1}{2}\sum_i w_i |b_i-C_n^b n_i|$
%%
% $L(A) = \sum_i w_i - tr\left(AB^T\right)$
%%
% where:
%%
% $B=\sum_i w_i b_i n_i^T$
%% 
% Attitude matrix that maximizes $tr(AB^T)$ minimizes the loss function $L(A)$
%%
% $tr\left(AB^T\right) = q^TKq$
%%
% where:
%%
% $K = \left[\begin{array}{cc}B+B^T-tr(B)I & \sum_i w_i b_i \times n_i \\
% \sum_i w_i \left(b_i \times n_i\right)^T & tr(B)\end{array}\right]$
%%
% q-method finds the optimal quaternion as the normalized eigenvector of K
% with the largest eigenvalue:
%%
% $Kq_{opt}=\lambda_{max} q_{opt}$

% Correct vector size
%fb = reshape(fb, [], 1); fn = reshape(fn, [], 1); 
%mb = reshape(mb, [], 1); mn = reshape(mn, [], 1);

%K matrix
%B = wf*fb*fn' + wm*mb*mn';
% Init la B matrix
[n,N] = size(fb);
B     = zeros(n,n);
% Calculate B
for i = 1:N
    B = B + w(i) * fb(:, i) * mn(:, i)';
end

K11 = B+B'-eye(3)*trace(B);
K22 = trace(B);

% Inicializar el vector z como un vector de ceros
K12 = zeros(3, 1);
% Calcular z de forma vectorial
for i = 1:N
    K12 = K12 + w(i) * cross(fb(:, i), mn(:, i));
end
%K12 = wf*cross(fb,fn)+wm*cross(mb,mn);
K21 = K12';
K = [K11, K12; K21, K22];

%Find eigenvalues and eigenvectors
[V,D] = eig(K);
lambdas = diag(D);

%Look for the largest eigenvalue
index = lambdas == max(lambdas);
q_opt = V(:,index);

%Transform q_opt to the Matlab's standard format
q = [q_opt(4), q_opt(1:3,1)'];
%Normalize quaternion
q = q/sqrt(q*q');

%As q-method yelds Cnb matrix
q = quatconj(q);
end