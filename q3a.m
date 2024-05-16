clear all
k=0;
t=linspace(0,2,50);
for i=-2:0.1:2
        x=(1-3*rho_0(i)^2).*t+i;
        plot(x,t);
        xlabel('x');
        ylabel('t')
        hold on
end

function rho=rho_0(x)
    if x<0
        rho=0.1;
    end
    if x>1
        rho=0.9;
    end
    if (x>=0 && x<=1)
        rho=0.1+0.8*x;
    end
end


