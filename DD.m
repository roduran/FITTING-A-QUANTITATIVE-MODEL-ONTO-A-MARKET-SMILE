function [ vol_DD, regression ] = DD( Spot, spotbarra, eta, omega  )
vol_DD=zeros(2663,1);


for i=1:2663
    vol_DD(i,1)=eta*(omega+(1-omega).*(spotbarra./Spot(i,1)));
end

x=1:2663;
x=x';
X=[x, ones(2663,1)];
coef=pinv(X)*vol_DD;
regression=polyval(coef,x);

end

