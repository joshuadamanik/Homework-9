function [cost,cnstr] = costt (p1)


   x1=p1(1);
   x2=p1(2);
   x3=p1(3);
   x4=p1(4);
   x5=p1(5);
   x6=p1(6);
   x7=p1(7);
   x8=p1(8);


   g1= -(1-0.0025*(x4+x6))*1000000;
   g2= -(1-0.0025*(x5+x7-x4))*1000000;
   g3= -(1-0.01*(x8-x5))*1000000;
   g4= -(x1*x6-833.33252*x4-100*x1+83333.333);  
   g5= -(x2*x7-1250*x5-x2*x4+1250*x4); 
   g6= -(x3*x8-1250000-x3*x5+2500*x5); 

%   g1= -(1-0.0025*(x4+x6));
%   g2= -(1-0.0025*(x5+x7-x4));
%   g3= -(1-0.01*(x8-x5));
%   g4= -(x1*x6-833.33252*x4-100*x1+83333.333)/1000000;
%   g5= -(x2*x7-1250*x5-x2*x4+1250*x4)/1000000;
%   g6= -(x3*x8-1250000-x3*x5+2500*x5)/1000000;

   cnstr=[g1 g2 g3 g4 g5 g6];

% cnstr = [ g1 g2 ... h1 h2 ....]   inequality constriants first 

   aa=x1+x2+x3;

   cost = aa;

