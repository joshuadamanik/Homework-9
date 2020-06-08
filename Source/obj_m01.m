function [cost,cnstr] = costfunc (p1)


   x1=p1(1);
   x2=p1(2);
   x3=p1(3);
   x4=p1(4);
   x5=p1(5);
   x6=p1(6);
   x7=p1(7);
   x8=p1(8);
   x9=p1(9);
   x10=p1(10);
   x11=p1(11);
   x12=p1(12);
   x13=p1(13);

   aa=5*(x1+x2+x3+x4)-5*(x1^2+x2^2+x3^2+x4^2);
   bb=-(x5+x6+x7+x8+x9+x10+x11+x12+x13);   

   cost = aa+bb;

  
   g1= (2*x1+2*x2+x10+x11-10);
   g2= (2*x1+2*x3+x10+x12-10);
   g3= (2*x2+2*x3+x11+x12-10);
   g4= (-8*x1+x10);
   g5= (-8*x2+x11);
   g6= (-8*x3+x12);
   g7= (-2*x4-x5+x10);
   g8= (-2*x6-x7+x11);
   g9= (-2*x8-x9+x12);

   cnstr=[g1 g2 g3 g4 g5 g6 g7 g8 g9];

% cnstr = [ g1 g2 ... h1 h2 ....]   inequality constriants first 