function tests = test_eta_c_inverse
% 运行方式：results=runtests('test_eta_c_inverse.m'); assertSuccess(results)
tests=functiontests(localfunctions);
end

function testSourceRegression(t)
% 在多个几何参数和正则化强度下，将独立 MATLAB 算子与回放结果
% 和只读参考夹具进行比较。
f=load(fullfile(fileparts(mfilename('fullpath')),'test_reference.mat'));
for k=1:numel(f.cases)
    ref=f.cases{k}; geo=ref.geometry;
    o=eta_c_inverse(geo(1),geo(2),geo(3),'Alpha',10.^f.logs);
    c=o.operator;
    ids=[1,round((numel(c.m2)+1)/2),numel(c.m2)];
    verifyEqual(t,c.m2(ids).',ref.sample_m2,'AbsTol',1e-11);
    verifyEqual(t,c.s,ref.s(:),'AbsTol',2e-13);
    verifyLessThan(t,norm(c.Khat(ids,:)-ref.Khat,'fro')/norm(ref.Khat,'fro'),2e-8);
    verifyLessThan(t,norm(c.g(ids)-ref.g(:))/norm(ref.g),3e-12);
    verifyLessThan(t,norm(c.H-ref.H,'fro')/norm(ref.H,'fro'),3e-12);
    for name={'one_delta','two_delta'}
        model=o.(name{1});
        if strcmp(name{1},'one_delta'), rr=ref.one; else, rr=ref.two; end
        r=eta_c_spectrum(model,10.^f.logs);
        for j=1:numel(f.logs)
            actual=[r.residues(j,:),r.rho(j,:)];
            expected=[rr.residues(j,:),rr.rho(j,:)];
            rel=norm(actual-expected)/max(norm(expected),eps);
            fprintf('source geo=%d %s loga=%g relative_coeff_error=%.3g\n',k,name{1},f.logs(j),rel);
            % 极小 alpha 的逆问题对独立组装 K 时的舍入误差较敏感。
            if f.logs(j)<-12, tol=2e-4; else, tol=5e-8; end
            verifyLessThan(t,rel,tol);
        end
        checkAugmented(t,c,model);
    end
end
end

function checkAugmented(t,c,model)
% 验证紧凑 GSVD 回放与直接增广最小二乘系统
% [A; sqrt(alpha)L] \ [b; 0] 一致。
f=model.factors;p=f.pab;ni=numel(c.s)-2;
A=sqrt(c.weights).*[exp(-model.masses.^2./c.m2)./c.m2,c.Khat(:,2:end-1)];
b=sqrt(c.weights).*(c.g-c.Khat(:,[1,end])*c.boundary)+A(:,p+1:end)*f.q;
L=[zeros(ni,p),chol(c.H(2:end-1,2:end-1))];
for a=10.^[-20,-12,-8,-4,0]
    r=eta_c_spectrum(model,a);
    direct=[A;sqrt(a)*L]\[b;zeros(ni,1)];
    err=norm(direct-r.shifted_coefficients)/max(norm(direct),eps);
    fprintf('aug p=%d alpha=%g relative_coeff_error=%.3g\n',p,a,err);
    verifyLessThan(t,err,3e-4);
    rr=norm(A*r.shifted_coefficients-b)/f.data_norm;
    verifyEqual(t,rr,r.relative_residual,'AbsTol',2e-12);
    verifyEqual(t,r.penalty^2,r.rho*c.H*r.rho.','RelTol',1e-10);
end
end

function testOffGridAndClosure(t)
% 检查精确端点处理、帽函数核闭合关系和最后的短单元。
f=load(fullfile(fileparts(mfilename('fullpath')),'test_reference.mat'),'off');
c=eta_c_operator(2.003,2.027,20.037); r=f.off;
verifyEqual(t,c.m2,r.m2(:),'AbsTol',1e-13);
verifyEqual(t,c.s,r.s(:),'AbsTol',1e-13);
verifyEqual(t,sum(c.weights),.024,'AbsTol',1e-14);
verifyLessThan(t,norm(c.g-r.g(:))/norm(r.g),3e-12);
verifyLessThan(t,norm(c.Khat-r.Khat,'fro')/norm(r.Khat,'fro'),2e-8);
verifyEqual(t,c.H,r.H,'AbsTol',1e-10);
verifyEqual(t,sum(c.Khat,2),exp(-c.s(1)./c.m2)-exp(-c.s(end)./c.m2),'AbsTol',2e-11);
verifyEqual(t,diff(c.s(1:end-1)),.1*ones(numel(c.s)-2,1),'AbsTol',1e-13);
verifyEqual(t,c.s(end),20.037);
end

function testSweepReplayAndNoSelection(t)
% MATLAB API 暴露完整 alpha 族，从不自动选择 alpha。
o=eta_c_inverse(2,2.1,20); % default 2001-alpha sweep, underdetermined data allowed
verifyEqual(t,numel(o.alpha),2001);
verifyEqual(t,o.alpha([1,end]),[1e-20,1]);
verifyEmpty(t,o.selected_alpha);
for name={'one_delta','two_delta'}
    m=o.(name{1}); r=eta_c_spectrum(m,o.alpha([1,1001,end]),[0,14.99,20,21]);
    verifyEqual(t,r.residues,m.residues([1,1001,end],:),'AbsTol',1e-10);
    verifyEqual(t,r.rho(:,1),zeros(3,1));
    verifyEqual(t,r.rho(:,end),repmat(o.operator.boundary(2),3,1));
    verifyEqual(t,r.rho_query(:,1:2),zeros(3,2));
    verifyTrue(t,all(r.rho_query(:,4)>0));
    verifyTrue(t,all(isnan(m.f(m.residues(:,1)<=0))));
end
% 单 alpha 输入：不应依赖极小值搜索、索引或自动选参假设。
a=eta_c_inverse(2,3,20,'Alpha',1e-4);
verifyEqual(t,numel(a.one_delta.f),1);
verifyEmpty(t,a.one_delta.selected_alpha);
end

function testInvalidInputs(t)
% 非法窗口、尺度和 alpha 网格必须显式报错。
bad={@()eta_c_inverse(NaN,80,30),@()eta_c_inverse(2,2,30), ...
    @()eta_c_inverse(1,80,30),@()eta_c_inverse(2,101,30), ...
    @()eta_c_inverse(2,80,14.99),@()eta_c_inverse(2,80,15.01), ...
    @()eta_c_inverse(2,80,30,'Alpha',[1,0]), ...
    @()eta_c_inverse(2,80,30,'Alpha',[1,1]), ...
    @()eta_c_inverse(2,80,30,'Alpha',[]), ...
    @()eta_c_inverse(2,80,30,'Alpha',Inf)};
for k=1:numel(bad)
    didThrow=false;
    try
        bad{k}();
    catch
        didThrow=true;
    end
    verifyTrue(t,didThrow,sprintf('invalid case %d must throw',k));
end
end
