q = 10;
l = 4;

n = 4;
n_el = zeros(1,n);
t = zeros(1,n);
M = @(q,l,x) -q*x.^2/2 + q*l/2*x - q*l^2/12;
v = @(q,l,x) -q*x.^4/24 + q*l/2*x.^3/6 - q*l^2/12*x.^2/2;
 
for j = 1:n
    tic
    dx = 10^(-j+1);
    x = 0:dx:l;
    fx = M(q,l,x);
    vx = v(q,l,x);
    for i = 1:2
        fx = integraleIndefinito(fx,x,0);
    end
    dif = (fx - vx)/max(abs(vx));
    n_el(j) = length(x)
    t(j) = toc
    figure(1)
    subplot(2,n,j)
    hold on
    plot(x,vx,'red')
    plot(x,fx)
    subplot(2,n,n+j)
    plot(x,dif)
    drawnow
end

