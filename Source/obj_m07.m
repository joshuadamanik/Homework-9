function [cost,cnstr] = costt (p1)


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

  
 
   g1= -(105-4*x1-5*x2+3*x7-9*x8);
   g2= -(-3*(x1-2)^2-4*(x2-3)^2-2*x3^2+7*x4+120);
   g3= -(-10*x1+8*x2+17*x7-2*x8);
   g4= -(-x1^2-2*(x2-2)^2+2*x1*x2-14*x5+6*x6);
   g5= -(8*x1-2*x2-5*x9+2*x10+12);
   g6= -(-5*x1^2-8*x2-(x3-6)^2+2*x4+40);
   g7= -(3*x1-6*x2-12*(x9-8)^2+7*x10);
   g8= -(-0.5*(x1-8)^2-2*(x2-4)^2-3*x5^2+x6+30);


   cnstr=[g1 g2 g3 g4 g5 g6 g7 g8];

   % cnstr = [ g1 g2 ... h1 h2 ....]   inequality constriants first 


   aa=x1^2+x2^2+x1*x2-14*x1-16*x2+(x3-10)^2+4*(x4-5)^2+(x5-3)^2;
   bb=2*(x6-1)^2+5*x7^2+7*(x8-11)^2+2*(x9-10)^2+(x10-7)^2+45;

   cost = aa+bb;
  
