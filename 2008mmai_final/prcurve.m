function prcurve()
p1 = [1/1,1/2,1/3,1/4,1/5];
r1 = [1/1,1/1,1/1,1/1,1/1];

p2 = [0/1,0/2,0/3,1/4,1/5];
r2 = [0/1,0/1,0/1,1/1,1/1];

p3 = [0/1,0/2,1/3,1/4,1/5];
r3 = [0/1,0/1,1/1,1/1,1/1];

p4 = [0/1,0/2,0/3,1/4,1/5];
r4 = [0/1,0/1,0/1,1/1,1/1];

p5 = [1/1,1/2,1/3,1/4,1/5];
r5 = [1/1,1/1,1/1,1/1,1/1];

plot(r1,p1,'o-',r2,p2,'x-',r3,p3,'s-',r4,p4,'.-',r5,p5,'*-');
legend('Cloth','Black Jacket','Scarf','Dress','White Jacket');
title('PR-Curve of different item query');
xlabel('Recall');
ylabel('Precision');
saveas(gcf,'prcurve.png','png');