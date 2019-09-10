%Matlab function for finding area from randomly distributed points

function area=poly_area(x,y)

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

    area=0;
    %finding area of the polygon
    for i=1:length(P)-1

        X1=P(1,i); X2=P(1,i+1);
        Y1=P(2,i); Y2=P(2,i+1);

        area=area+(X1*Y2-X2*Y1);

    end


    area=0.5*(area);
    
end