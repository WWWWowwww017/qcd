function out = eta_c_inverse(M2min, M2max, Lambda, varargin)
%ETA_C_INVERSE 独立的 A-A n=0 LO+G2+G3 逆问题，重整化尺度 mu=2 GeV。
% out = eta_c_inverse(M2min,M2max,Lambda) 返回单双极点两种模型。
% 可选参数 Alpha 为正且严格递增的向量，默认值为 10.^(-20:.01:0)。
% 函数不自动选择 alpha；eta_c_spectrum(out.one_delta,alpha) 可回放任意 alpha。
% 归一化、来源和适用边界见 README.md。

% 三个输入定义 Borel 窗口和连续谱分离尺度。函数有意返回完整 alpha
% 族，而不是静默地选择某一个正则化强度。
validateattributes(M2min,{'numeric'},{'real','finite','scalar','>=',1.5},mfilename,'M2min');
validateattributes(M2max,{'numeric'},{'real','finite','scalar','<=',100},mfilename,'M2max');
validateattributes(Lambda,{'numeric'},{'real','finite','scalar','>',14.99},mfilename,'Lambda');
if M2max <= M2min
    error('eta_c:Window','Require 1.5 <= M2min < M2max <= 100 (GeV^2).');
end
p = inputParser;
p.addParameter('Alpha',10.^linspace(-20,0,2001),@validAlpha);
p.parse(varargin{:});
alpha = double(p.Results.Alpha(:).');

% eta_c_operator 组装 Borel 数据 g(M^2)、帽函数基核、H1 惩罚矩阵，
% 以及 rho(s) 的固定边界值。
c = eta_c_operator(double(M2min),double(M2max),double(Lambda));
if numel(c.s) < 3
    error('eta_c:Grid','Lambda must leave at least one interior hat node (Lambda > 15.09).');
end
out = struct('geometry',[M2min,M2max,Lambda], 'mu',2, ...
    'mc',1.1334710909,'G2',0.038,'G3',0.013,'s_h',14.99, ...
    'ds',0.1,'dM2',0.01,'alpha',alpha,'operator',c, ...
    'selected_alpha',[],'selection','none: caller chooses alpha');

% 两种模型使用相同的连续谱离散化。单极点模型保留一个不受惩罚的留数，
% 双极点模型再增加第二个留数。
out.one_delta = solveModel(c,2.9841,alpha);
out.two_delta = solveModel(c,[2.9841,3.51067],alpha);
end

function ok = validAlpha(a)
ok = isnumeric(a) && isreal(a) && isvector(a) && ~isempty(a) ...
    && all(isfinite(a(:))) && all(a(:)>0) && all(diff(a(:))>0);
end

function model = solveModel(c,masses,alpha)
pab = numel(masses);
ni = numel(c.s)-2;

% 消去两个固定的连续谱端点。q 是 H1 交叉项引起的平移，
% 因而自由内部变量是在给定边界附近求解的。
Hii = c.H(2:end-1,2:end-1);
q = Hii \ (c.H(2:end-1,[1,end])*c.boundary);
L = [zeros(ni,pab),chol(Hii)]; % MATLAB upper Cholesky: L'*L = Hii

% A 包含极点列和内部帽函数列。Borel 求积权重被吸收到 A 与 b 中，
% 从而数据范数采用欧氏范数。
poles = exp(-masses.^2 ./ c.m2) ./ c.m2;
A = sqrt(c.weights) .* [poles,c.Khat(:,2:end-1)];
if rank(A(:,1:pab)) < pab
    error('eta_c:PoleRank','Borel window cannot resolve the unpenalized pole columns.');
end
b = sqrt(c.weights).*(c.g-c.Khat(:,[1,end])*c.boundary) ...
    + A(:,pab+1:end)*q;

% 紧凑 GSVD 前先用 QR 去除高矩阵数据的零空间残差。
% 广义奇异值 (cc, mu) 将 Tikhonov 问题对角化。
[Q,R] = qr(A,0);
bq = Q'*b;
null2 = norm(b-Q*bq)^2;
[U,~,X,C,S] = gsvd(R,L,0);
cc = sqrt(sum(C.^2,1)).';
mu = sqrt(sum(S.^2,1)).';
% C'*U'*bq 等于 cc.*beta，其中也包括数据奇异值为零的情形。
numerator = C'*(U'*bq);
factors = struct('cc',cc,'mu',mu,'numerator',numerator, ...
    'invXt',X'\eye(size(X,1)),'q',q,'pab',pab, ...
    'boundary',c.boundary,'H',c.H,'s',c.s,'masses',masses, ...
    'R',R,'bq',bq,'null2',null2,'data_norm',norm(sqrt(c.weights).*c.g));
if any(cc.^2+mu.^2 == 0) || factors.data_norm == 0
    error('eta_c:Degenerate','Degenerate GSVD or zero OPE norm.');
end
model = struct('masses',masses,'factors',factors,'selected_alpha',[]);

% 分块计算标量扫描诊断量。完整谱函数由 eta_c_spectrum 按需重构，
% 以保持归档紧凑。
n = numel(alpha);
model.alpha = alpha(:);
model.residues = zeros(n,pab);
model.f = nan(n,1);
model.df_dlog10alpha = nan(n,1);
model.relative_residual = zeros(n,1);
model.penalty = zeros(n,1);
model.minrho = zeros(n,1);
model.negative_count = zeros(n,1);
for first = 1:128:n
    ids = first:min(first+127,n);
    r = eta_c_spectrum(model,alpha(ids));
    model.residues(ids,:) = r.residues;
    model.f(ids) = r.f;
    model.df_dlog10alpha(ids) = r.df_dlog10alpha;
    model.relative_residual(ids) = r.relative_residual;
    model.penalty(ids) = r.penalty;
    model.minrho(ids) = min(r.rho,[],2);
    model.negative_count(ids) = sum(r.rho<0,2);
end
end
