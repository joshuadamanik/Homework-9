function [cost,cnstr] = costt (p1)


   x1=p1(1);
   x2=p1(2);
   x3=p1(3);
   x4=p1(4);
   x5=p1(5);
   x6=p1(6);
   x7=p1(7);
  
 
   g1= -(127-2*x1^2-3*x2^4-x3-4*x4^2-5*x5);
   g2= -(282-7*x1-3*x2-10*x3^2-x4+x5);
   g3= -(196-23*x1-x2^2-6*x6^2+8*x7);
   g4= -(-4*x1^2-x2^2+3*x1*x2-2*x3^2-5*x6+11*x7);

   cnstr=[g1 g2 g3 g4];

% cnstr = [ g1 g2 ... h1 h2 ....]   inequality constriants first 

   aa=(x1-10)^2+5*(x2-12)^2+x3^4+3*(x4-11)^2;
   bb=10*x5^6+7*x6^2+x7^4-4*x6*x7-10*x6-8*x7;

   cost = aa+bb;

