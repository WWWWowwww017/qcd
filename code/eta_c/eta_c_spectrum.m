function r = eta_c_spectrum(model,alpha,squery)
%ETA_C_SPECTRUM 根据保存的因子回放选定的正 alpha 值。
% r=eta_c_spectrum(out.one_delta,alpha) 或 eta_c_spectrum(out.two_delta,alpha)。
% 可选参数 squery 在 [s_h,Lambda] 上插值连续谱，s_h 以下置零，
% Lambda 以上使用固定的 LO 微扰连续谱。Delta 极点单独返回质量和留数，
% 从不使用有限宽度峰近似。rho 的行对应 alpha，列对应 r.s。

% GSVD 解在紧凑模态中是对角的。每个模态乘以滤波因子
% cc/(cc^2 + alpha*mu^2)；未正则化的极点模态予以保留。
validateattributes(alpha,{'numeric'},{'real','finite','vector','positive','nonempty'});
if ~isstruct(model) || ~isscalar(model) || ~isfield(model,'factors')
    error('eta_c:Model','Pass one_delta or two_delta returned by eta_c_inverse.');
end
f=model.factors; alpha=double(alpha(:).');
den=f.cc.^2+f.mu.^2.*alpha;
z=f.invXt*(f.numerator./den);
dz=f.invXt*(-log(10)*f.numerator.*(f.mu.^2.*alpha)./den.^2);

% 将平移后的未知量还原为物理极点留数和连续谱节点，
% 撤销 q 后恢复固定端点值。
r=struct('alpha',alpha(:),'s',f.s,'masses',f.masses, ...
    'residues',z(1:f.pab,:).','shifted_coefficients',z, ...
    'rho',[repmat(f.boundary(1),1,numel(alpha)); ...
        z(f.pab+1:end,:)-f.q;repmat(f.boundary(2),1,numel(alpha))].');
r.f=nan(numel(alpha),1);
good=r.residues(:,1)>0;
r.f(good)=sqrt(r.residues(good,1));
r.df_dlog10alpha=nan(size(r.f));
r.df_dlog10alpha(good)=dz(1,good).'./(2*r.f(good));
r.relative_residual=sqrt(sum((f.R*z-f.bq).^2,1)+f.null2).'/f.data_norm;
r.penalty=sqrt(max(sum(r.rho.'.*(f.H*r.rho.'),1),0)).';
if nargin>=3
    validateattributes(squery,{'numeric'},{'real','finite','vector','nonempty','nonnegative'});
    sq=double(squery(:));
    rq=zeros(numel(alpha),numel(sq));
    % 只插值重构得到的低能连续谱。Lambda 以上使用固定的 LO 微扰尾部，
    % 阈值以下置零。
    mid=sq>=f.s(1) & sq<=f.s(end);
    rq(:,mid)=interp1(f.s,r.rho.',sq(mid),'linear').';
    high=sq>f.s(end); v=sqrt(1-4*1.1334710909^2./sq(high));
    rq(:,high)=repmat((3*v.*(1-v.^2/3)/(8*pi^2)).',numel(alpha),1);
    r.squery=sq; r.rho_query=rq;
end
end
