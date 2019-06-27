function [S] = SpaceNodes( Smax, Smin, N )
S=zeros(1,N);
ds=(Smax-Smin)/(N-1);
S(1,1)=Smin;
for i=1:N-1
    S(1,i+1) = S(1,i)+ds;
end

end

