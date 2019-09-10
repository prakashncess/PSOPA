%Matlab function for area of any random polynomial

function poly_plot(xx)

    %x is any points in x coordinates; 
    %y is any points in y coordinates;
    sz=length(xx);
    x=xx(1,1:sz/2);
    y=xx(1,(sz/2)+1:sz);
    P=[x;y];
    c=mean(P,2);
    d=P-c;
    th=atan2(d(2,:),d(1,:));
    [th, idx]=sort(th);
    P=P(:,idx);
    P=[P P(:,1)];
    fill(P(1,:),P(2,:),'.-r')
    hold on
    plot(x,y,'b*');
    set(gca,'YDir','reverse')
end