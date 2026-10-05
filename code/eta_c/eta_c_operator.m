function c = eta_c_operator(M2min,M2max,Lambda)
%ETA_C_OPERATOR 纯 MATLAB 固定 OPE、物理 s 网格帽函数及 H1 矩阵。
% 主要用于复现与测试，几何参数域与逆问题函数相同。

% m2 是 Borel 网格，s 是物理谱网格。两者都包含用户给定的精确端点；
% 最后一个谱单元可能短于名义 Delta-s 步长。
validateattributes(M2min,{'numeric'},{'real','finite','scalar','>=',1.5});
validateattributes(M2max,{'numeric'},{'real','finite','scalar','<=',100,'>',M2min});
validateattributes(Lambda,{'numeric'},{'real','finite','scalar','>',14.99});
M2min=double(M2min); M2max=double(M2max); Lambda=double(Lambda);
lattice = (1.5:0.01:100).';
m2 = [M2min;lattice(lattice>M2min+1e-10 & lattice<M2max-1e-10);M2max];
n = floor((Lambda-14.99)/0.1+1e-10);
s = 14.99+0.1*(0:n).';
if Lambda-s(end)>1e-10
    s(end+1)=Lambda;
else
    s(end)=Lambda;
end
if numel(s)<2 || any(diff(s)<=0)
    error('eta_c:Grid','Lambda too close to threshold for a resolvable cell.');
end
K = zeros(numel(m2),numel(s));
H = zeros(numel(s),numel(s));

% 对每个分段线性帽函数基，解析积分其与 exp(-s/M^2)/M^2 的乘积。
% 当 h/M^2 较小时使用级数展开，以避免消去误差。
for j=1:numel(s)-1
    h=s(j+1)-s(j); z=h./m2;
    small=abs(z)<1e-4;
    asc=zeros(size(z)); des=asc;
    t=z(small);
    asc(small)=0.5-t/3+t.^2/8-t.^3/30+t.^4/144;
    des(small)=0.5-t/6+t.^2/24-t.^3/120+t.^4/720;
    t=z(~small); e=exp(-t);
    asc(~small)=(1-(1+t).*e)./t.^2;
    des(~small)=(t-1+e)./t.^2;
    fac=exp(-s(j)./m2).*z;
    K(:,j)=K(:,j)+fac.*des;
    K(:,j+1)=K(:,j+1)+fac.*asc;
    % 在当前单元上组装有限元 H1 质量项与刚度项矩阵。
    H(j:j+1,j:j+1)=H(j:j+1,j:j+1) ...
        +[h/3+1/h,h/6-1/h;h/6-1/h,h/3+1/h];
end

% 两个端点值固定：强子阈值处为零，Lambda 处为 LO 微扰连续谱。
% OPE 给出数据向量。
[pert,g2,g3]=ope(Lambda,m2);
c=struct('m2',m2,'s',s,'weights',trapweights(m2), ...
    'Khat',K,'H',H,'boundary',[0;rhoPert(Lambda)], ...
    'g',pert+g2+g3,'pert',pert,'g2',g2,'g3',g3);
end

function r=rhoPert(s)
mc=1.1334710909;
r=zeros(size(s)); mask=s>4*mc^2;
v=sqrt(1-4*mc^2./s(mask));
r(mask)=3*v.*(1-v.^2/3)/(8*pi^2);
end

function [pert,g2,g3]=ope(Lambda,m2)
mc=1.1334710909; q=4*mc^2;

% 微扰项在速度变量中积分；D4/D6 凝聚项使用费曼参数 x。
% 这些是 README.md 记录的固定 A-A n=0 表达式。
v=linspace(0,sqrt(max(1-q/Lambda,0)),16001).';
sg=q./max(1-v.^2,q/Lambda); sg(end)=Lambda;
rw=rhoPert(sg).*trapweights(sg);
pert=zeros(size(m2)); g2=pert; g3=pert;
x=linspace(1e-8,1-1e-8,4001); xb=1-x;
wx=trapweights(x.');
a3=-(45/8)*x.*xb.*(x.^3+xb.^3)-(69/72)*x.^2.*xb.^2;
b3=-(3/4)*mc^2*x.*xb.*(x.^4+xb.^4) ...
    +mc^2*x.^2.*xb.^2.*(-(23/12)*(x.^2+xb.^2)+x.*xb/3+2);
c3=(2/5)*mc^4*x.*xb.*(x.^5+xb.^5);
for a=1:128:numel(m2)
    ids=a:min(a+127,numel(m2)); mm=m2(ids);
    % 分块处理 Borel 点，使临时指数数组保持在可控大小。
    pert(ids)=(exp(-sg.'./mm)./mm)*rw;
    E=exp(-mc^2./(x.*xb.*mm));
    i2=E.*(1./(2*mm.^2)-mc^2*(x.^3+xb.^3)./(2*mm.^3.*x.^2.*xb.^2));
    i3=E.*(a3./(2*mm.^3.*x.^3.*xb.^3) ...
        +b3./(6*mm.^4.*x.^4.*xb.^4)+c3./(24*mm.^5.*x.^5.*xb.^5));
    g2(ids)=0.038*(i2*wx)/(6*pi);
    g3(ids)=0.013*(i3*wx)/(4*pi)^2;
end
end

function w=trapweights(x)
h=diff(x);
w=[h(1)/2;(h(1:end-1)+h(2:end))/2;h(end)/2];
end
