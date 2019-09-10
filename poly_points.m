%Matlab function for arranging points in counterclockwise from randomly
%distributed points

function [x1,y1]=poly_points(x,y)

    %x is any points in x coordinates; 
    %y is any points in y coordinates;

    P=[x;y];
    c=mean(P,2);
    d=P-c;
    th=atan2(d(2,:),d(1,:));
    [th, idx]=sort(th);
    P=P(:,idx);
    P=[P P(:,1)];
    %fill(P(1,:),P(2,:),'.-r')
    %hold on
    %plot(x,y,'b*')
    
    %x1 and y1 are arraged anticlockwise points for points(x,y)
    x1=P(1,1:end-1);
    y1=P(2,1:end-1);
    
end