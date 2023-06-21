function FitEllipse2Cloud(p)

if nargin == 0
    p = randn(1000,2);
    p = p.*[1 5];
    p = p*[cos(pi/4) sin(pi/4);-sin(pi/4) cos(pi/4)];
end

mu = mean(p,1);
p  = p-mu;

% C = (p'*p)/size(p,1)
C = cov(p);

[V,D] = eigs(C);

Maj = sqrt(D(2,2));
Min = sqrt(D(1,1));

alpha = atan2d(V(1,2),V(1,1));

p = p+mu;
plot(p(:,1),p(:,2),'b.')
axis image

hold on
% alpha;
drawellipse('SemiAxes',[2*Maj 2*Min],'RotationAngle',-alpha,'Center',mu,...
    'markersize',eps)

