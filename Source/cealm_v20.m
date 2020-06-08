%---------------------------------------------------------------------------------
%  CEALM (Co-Evolution Augmented Lagrangian Method)  Version 2.0
%---------------------------------------------------------------------------------
%
%%%  Ver 2.0
%%%  - Equality constraints and inequality constraints can be handled simultaneously
%%%  - Initial bounds of multipliers are set automatically 
%%%
%%%  Ver 1.2 
%%%%  - Man-to-man strategy calculation modified (fewer individual ranked)
%%%
%%%  Ver 1.1 
%%%  - Print option modified
%%%  - Man-to-man strategy included
%%%
%
%  Function name : cealm_v20.m  
%
%  by TAHK Min-Jea (mjtahk@fdcl.kaist.ac.kr)
%
%  Ver 2.0 - September 9, 2004
%  Ver 1.0 - May 29, 1998  
%  
%
%
%  Objective:
%
%       To solve constrained optimization problems using co-evoultuion 
%
%  Algorithm:
%
%	1. Co-evolution algorithms for saddle-point problems  
%	2. Equality and Inequality Constraints 
%	3. ES for recombination & mutation
%
%
%  Variables:  (U=user defined)
%
%  U    NumParentX(Y)    : number of parent population   
%  U    NumOffspringX(Y) : number of offspring population   
%  U    NumParamX        : number of parameters
%       NumParamY        : number of multipliers
%  U    NumIneq          : number of inequality constraints
%  U    NumEq            : number of equality constraints
%       ParentPopX(Y)    : parameter values of parent population
%       ParentSigmaX(Y)  : standard deviations of parent population
%       OffspringPopX(Y) : paramter values of offspring population
%       SigmaX(Y)        : standard deviations of offspring population
%  U    ParamMaxX        : upper bounds of search space for parameters
%  U    ParamMinX        : lower bounds of search space for parameters
%       ParamMaxY        : upper bounds of search space for multipliers
%       ParamMinY        : lower bounds of search space for multipliers
%  U    Tolerance        : tolerance of constraint violation     
%       CV               : inequality constraint violation     
%  U    varX(Y)          : initial standard deviation relative to search interval   
%       Gencount         : generation count
%  U    MaxGen       	 : maximum number of generation
%  U    SigInitX(Y)      : initial sigma bound (during InitFlag==0)
%  U    SigMinX(Y)       : final sigma bound (during InitFlag==1)
%
%  Flags:
%       InitFlag    = 1 if all the constraints are satisfied for the first time by the best X  
%       CnstrFlag   = 1 if all the constraints are satisfied
%       ImprvFlag    = 1 if the cost is improved while satisfying constraints
%
%  istrategy = 0   both players use security strategy 
%  istrategy = 1   the leader(the multipliers) uses security strategy but the follower uses man-to-man strategy
%
%  Screen Outputs:
%     Run = run count for Monte-Carlo sim
%     G   = generation number
%     BestF = best cost obtained during the evolution
%     BestV = constraint violation for the case of BestF
%     X, Y  = parameter and multipler values for the case of BestF
%     CostF = cost value for the best individual X at the final generation
%     Vsum = constraint violation of the best individual at the final generation
%     AL_X, AL_Y = ALM values for the best individuals X and Y at the final generation
%     X, Y = values of the best individuals X and Y at the final generation
 
%--------- variance parameters ----------------------------------------------

   varX		 = 0.2;                 % initial sigma devided by interval    
   varY      = 0.2;

   KK        = 1;

      
%---------- initial bounds of Lagrange multipliers

   NumParamY = NumIneq + NumEq;
   
   if(NumIneq ~= 0)
       for k=1:NumIneq
          ParamMaxY(k) = 10;
          ParamMinY(k) = 0;     
       end
   end    
   if(NumEq ~= 0)
       for k=1:NumEq
          ParamMaxY(NumIneq+k) = 10;
          ParamMinY(NumIneq+k) = -10;     
       end
   end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Gencount 	 = 1 ;  	                  % Generation counter

   InitFlag      = 0 ;

   IntervalX     = ParamMaxX - ParamMinX;
   IntervalY     = ParamMaxY - ParamMinY;

   NumEval       = NumOffspringX/10^6;
 
   TtlParam      = NumParamX + NumParamY;

   taupX         = KK/(sqrt(2*TtlParam));	   % tau and tau prime according to
   tauX          = KK/(sqrt(2*sqrt(TtlParam)));   % Ba"ck's book

   taupY         = KK/(sqrt(2*TtlParam));	
   tauY          = KK/(sqrt(2*sqrt(TtlParam)));	

   ParentPopX    = zeros(NumParentX,NumParamX);
   ParentPopY    = zeros(NumParentY,NumParamY);
   ParentSigmaX  = zeros(NumParentX,NumParamX);
   ParentSigmaY  = zeros(NumParentY,NumParamY);

   Dummy         = ones(1,NumParamY);

   Tol           = Tolerance*ones(1,NumParamY);              

   PrevF          = 1000000000000000.;
   PrevV          = 1000000000000000.;

%---------------------------------------------------
% Initial Population: 
%    Uniform distrubution within the search space
%---------------------------------------------------
  
  rand('seed',sum(100*clock));                % Used 'seed' instead of 'state' for 
  randn('seed',sum(100*clock));               % compatability of Matcom3 to Matlab4

%   VarIntervalX  = varX*IntervalX;
%   VarIntervalY  = varY*IntervalY;
  

  for i=1:NumParamX 
        TempX(i)  = 0.2D0;
        TempXf(i)  = 0.1D-12;  
        VarIntervalX(i) = TempX(i)*IntervalX(i);
  end
    
  for i=1:NumParamY
        TempY(i)  = 0.2D0;
        TempYf(i)  = 0.1D-12;
        VarIntervalY(i) = TempY(i)*IntervalY(i);
  end

  for i=1:NumParentX  
     for j=1:NumParamX
         ParentPopX(i,j) = rand*(ParamMaxX(j)-ParamMinX(j)) + ParamMinX(j);
         ParentSigmaX(i,j)= VarIntervalX(j);
     end
  end

  for i=1:NumParentY   
     for j=1:NumParamY
        ParentPopY(i,j) = rand*(ParamMaxY(j)-ParamMinY(j)) + ParamMinY(j);
        ParentSigmaY(i,j)= VarIntervalY(j);
     end
  end   

SigInitX      = 0.000001*ones(1,NumParamX);
SigInitY      = 0.000001*ones(1,NumParamY);
  
%  SigInitX = VarIntervalX;
%  SigInitY = VarIntervalY;

SigMinX       = 0.000001*ones(1,NumParamX);
SigMinY       = 0.000001*ones(1,NumParamY);

%--------- mutation lower bounds ---------------------------------------------

   alpha  = 10^(-1/DeciGen);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% >>======= Iteration Starts Here =======<< %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPrint = 1;
while (Gencount <= MaxGen)  

%-----------------------------------------------------
% Recombination Process  (evolutionary strategy)
%-----------------------------------------------------


 for i=1:NumOffspringX
    SX = ceil(rand*NumParentX); % Random Number Generation
    TX = ceil(rand*NumParentX); % for Recombination
    for j=1:NumParamX
       OffspringPopX(i,j) = ParentPopX(SX,j);
       SigmaX(i,j)        = (ParentSigmaX(SX,j)+ParentSigmaX(TX,j))/2;
    end
 end

 for i=1:NumOffspringY
    SY = ceil(rand*NumParentY);
    TY = ceil(rand*NumParentY);
    for j=1:NumParamY
       OffspringPopY(i,j) = ParentPopY(SY,j);
       SigmaY(i,j)        = (ParentSigmaY(SY,j)+ParentSigmaY(TY,j))/2;
    end
 end

 
%-------------------------------------------------------------------------
% Mutation Process  (evolutionary strategy with an annealing feature)
%-------------------------------------------------------------------------

for i=1:NumParamX
      if TempX(i) > TempXf(i)
	     TempX(i) = alpha*TempX(i);
      else
	     TempX(i) = TempXf(i);
      end
end

for  i=1:NumParamY
      if TempY(i) >  TempYf(i) 
	     TempY(i) = alpha*TempY(i);
      else
	     TempY(i) = TempYf(i);
      end
end

  
for i=1:NumParamX
	  SigmaMinX(i) = TempX(i)*IntervalX(i);
end
for i=1:NumParamY
	  SigmaMinY(i) = TempY(i)*IntervalY(i);
end

%if InitFlag == 0
%  SigBoundX = SigInitX; 
%  SigBoundY = SigInitY;   
%else
%  SigBoundX = SigMinX;
%  SigBoundY = SigMinY;
%end



for i=1:NumOffspringX
   randnX = randn;
   for j=1:NumParamX
	SigmaX(i,j) = SigmaX(i,j)*(exp(taupX*randnX +tauX*randn));
        if SigmaX(i,j) < SigmaMinX(j)
            SigmaX(i,j) = SigmaMinX(j);
        end
        if SigmaX(i,j) > VarIntervalX(j)
            SigmaX(i,j) = VarIntervalX(j);
	end
	OffspringPopX(i,j) = OffspringPopX(i,j) + SigmaX(i,j)*randn;
   end
end

for i=1:NumOffspringY
   randnY = randn;
   for j=1:NumParamY
	SigmaY(i,j) = SigmaY(i,j)*(exp(taupY*randnY + tauY*randn));
        if SigmaY(i,j) < SigmaMinY(j)
            SigmaY(i,j) = SigmaMinY(j);
        end
        if SigmaY(i,j) > VarIntervalY(j)
            SigmaY(i,j) = VarIntervalY(j);
	end
	OffspringPopY(i,j) = OffspringPopY(i,j) + SigmaY(i,j)*randn;
   end
end


%--------------------------------------------
%  search space restriction (softwalls)
%--------------------------------------------

for i=1:NumOffspringX
   for j=1:NumParamX
      violate1 = OffspringPopX(i,j)-ParamMaxX(j);
      violate2 = ParamMinX(j)-OffspringPopX(i,j);
      if violate1 >0 
         OffspringPopX(i,j) = ParamMaxX(j) - rand*min(violate1,IntervalX(j)) ;
      end
      if violate2 >0 
         OffspringPopX(i,j) = ParamMinX(j) + rand*min(violate2,IntervalX(j)) ;
      end
   end
end

%-------------- lower bound check only for inequality constraints -------------- 

for i=1:NumOffspringY
   for j=1:NumIneq                        
      violate2 = ParamMinY(j)-OffspringPopY(i,j);
      if violate2 >0 
         OffspringPopY(i,j) = ParamMinY(j)+rand*min(violate2,IntervalY(j)) ;
      end
   end
end



%--------------------------------------------
% Matching Process (full match) 
%--------------------------------------------


for i=1:NumOffspringX

  [costf,cnstr] = feval(CostDef,OffspringPopX(i,:));

  for j=1:NumOffspringY
    AL1 = 0;
    AL2 = 0;
    AL3 = 0;
    if(NumIneq ~= 0)
        for k=1:NumIneq
           AL1 = AL1 + (max(cnstr(k)+0.5*OffspringPopY(j,k)/rho, 0))^2;
           AL2 = AL2 + OffspringPopY(j,k)^2;
        end
        AL1 = rho*AL1;
        AL2 = - 0.25*AL2/rho; 
    end
    if(NumEq ~= 0)
        for k=1:NumEq
           AL3 = AL3 + OffspringPopY(j,NumIneq+k)*cnstr(NumIneq+k) + rho*cnstr(NumIneq+k)^2;  
        end
    end
    ALV = costf + AL1 + AL2 + AL3;
    CostOffspringX(i,j) = ALV;
    CostOffspringY(j,i) = ALV; 
  end
end


%--------------------------------------------------------
%  Fitness evaluation and ordering
%--------------------------------------------------------


if(istrategy==0)    % Y & X = security
    [maxcost,mxcind] = max(CostOffspringX');  
    [mincost,mncind] = min(CostOffspringY');
    [SortedCostX,IndexX] = sort(maxcost);    % smallest first
    [SortedCostY,IndexY] = sort(-mincost);
end

if(istrategy==1)    % Y=security, X=man-to-man
    
     [mincost,mncind] = min(CostOffspringY');
     [SortedCostY,IndexY] = sort(-mincost);    % largest first

     xlarge=10^100; 
    for j=1:NumParentX
 %    for j=1:NumOffspringY
        ind = IndexY(j);
        yrow = CostOffspringY(ind,1:NumOffspringX);
        [minrow,mnix]=min(yrow);
        IndexX(j)=mnix;
        SortedCostX(j)=minrow;
        CostOffspringY(:,mnix)=xlarge;
     end 
  end
  
  SortedCostX = SortedCostX';
  SortedCostY = -SortedCostY';

  
  SortedPopulationX 	 = OffspringPopX(IndexX,:);
  SortedSigmaX		 = SigmaX(IndexX,:);

  SortedPopulationY 	 = OffspringPopY(IndexY,:);
  SortedSigmaY		 = SigmaY(IndexY,:);


%------------------------------------------------------
% Best individuals of the current generation 
%------------------------------------------------------
    
   [BestCostX,BestIndividualX] = min(SortedCostX);		
   BestParamX = SortedPopulationX(BestIndividualX,:);
%   WorstParamX = SortedPopulationX(NumOffspringX,:);
   BestSigmaX = SortedSigmaX(BestIndividualX,:);
%   WorstSigmaX = SortedSigmaX(NumOffspringX,:);

   [BestCostY,BestIndividualY] = max(SortedCostY);		
   BestParamY = SortedPopulationY(BestIndividualY,:);
%   WorstParamY = SortedPopulationY(NumOffspringY,:);
   BestSigmaY = SortedSigmaY(BestIndividualY,:);
%   WorstSigmaY = SortedSigmaY(NumOffspringY,:);


%--------------------------------------------------
% Output to Display & File
%--------------------------------------------------

%  TtlEval = NumEval*Gencount;

  [CostF,Cnstr] = feval(CostDef,BestParamX);
  
  if(NumEq ~= 0) 
     for k=1:NumEq
        Cnstr(NumIneq+k)=abs(Cnstr(NumIneq+k));
     end
  end    
  Vx  = max(Cnstr,0);
  Vsum = sum(Vx);
  
  C1 = Cnstr - Tol;
  V1   = max(C1,0);
  CV  = sum(V1);


  if CV <= 0
    CnstrFlag = 1;
  else
    CnstrFlag = 0;
  end

  if InitFlag == 0  & CnstrFlag ==1              % set BestF if constraint is satisfied 
    InitFlag = 1;                                % for the first time 
    PrevF = CostF;
    PrevV = Vsum;
    PrevPX = BestParamX;
    PrevPY = BestParamY;
  end

  if InitFlag ==1 & CnstrFlag ==1  & CostF < PrevF      % check if cost improved   
    ImprvFlag = 1;
  else
    ImprvFlag=0;
  end

  if ImprvFlag == 1
     PrevF = CostF;
     PrevV = Vsum;
     PrevPX = BestParamX;
     PrevPY = BestParamY;
  end

%------  Write to the output file if improved ------------------

  if ImprvFlag == 1 || (InitFlag == 0  & CnstrFlag ==1 ) 
    fprintf(outfile, '\n Run=%3d   G=%3d', irun, Gencount) ;
    fprintf(outfile, '\n CostF =%10.7f', CostF) ;
    fprintf(outfile,'\n');
    if NumIneq ~= 0 
        for i=1:NumIneq
           fprintf(outfile,' G%d=%10.7f',i, Cnstr(i));
        end
        fprintf(outfile,'\n');
    end
    if NumEq ~= 0 
       for i=1:NumEq
           fprintf(outfile,' H%d=%10.7f',i, Cnstr(i+NumIneq));
       end
       fprintf(outfile,'\n');
    end
    for i=1:NumParamX
      fprintf(outfile,' X%d=%10.7f',i, BestParamX(i));
    end
    fprintf(outfile,'\n');
    for i=1:NumParamY
      fprintf(outfile,' Y%d=%10.7f',i, BestParamY(i));
    end
    fprintf(outfile,'\n');
  end   
     
%-------- Print at every nPrint generation -------------------------------------------

   if iPrint ==  nPrint
     fprintf('\n Run=%3d G=%3d ', irun, Gencount);
     fprintf('\n CV    =%10.7f',CV);
     for i=1:NumParamY
       fprintf('   G%d=%10.7f',i, Cnstr(i));
     end
     fprintf('\n CostX =%10.7f', BestCostX) ;
     for i=1:NumParamX
       fprintf('   X%d=%10.7f',i, BestParamX(i));
     end
     fprintf('\n CostY =%10.7f', BestCostY) ;
     for i=1:NumParamY
       fprintf('   Y%d=%10.7f',i, BestParamY(i));
     end
     fprintf('\n');
     iPrint = 0;
   end
   iPrint = iPrint + 1;
   
 
%--------Screen output at the end of evolution --------------

  if(Gencount >= MaxGen)
     fprintf('\n Evolution Summary: Run=%3d G=%3d ', irun, Gencount);

    if(InitFlag == 1)
%       fprintf('\n BestF=%10.7f BestV=%10.7f', PrevF, PrevV);
%       fprintf('\n');
%       for i=1:NumParamX
%          fprintf('   X%d=%10.7f',i, PrevPX(i));
%        end
%        fprintf('\n');
%        for i=1:NumParamY
%          fprintf('   Y%d=%10.7f',i, PrevPY(i));
%        end
     else
        fprintf('   Unable to find a feasible solution!');
     end
    
     fprintf('\n CostF=%10.7f Vsum=%10.7f AL_X=%10.7f AL_Y=%10.7f', CostF, Vsum, BestCostX, BestCostY);
     fprintf('\n');
     for i=1:NumParamX
       fprintf('   X%d=%10.7f',i, BestParamX(i));
     end
     fprintf('\n');
     for i=1:NumParamY
       fprintf('   Y%d=%10.7f',i, BestParamY(i));
     end
     fprintf('\n');
      
 end


%--------------------------------------
% New Parent for next generation 
%--------------------------------------

 ParentPopX(1:NumParentX,:)   = SortedPopulationX(1:NumParentX,:);
 ParentPopY(1:NumParentY,:)   = SortedPopulationY(1:NumParentY,:);
 ParentSigmaX(1:NumParentX,:) = SortedSigmaX(1:NumParentX,:);
 ParentSigmaY(1:NumParentY,:) = SortedSigmaY(1:NumParentY,:);

%  check_transient;
 graph_data;
%  if done == true, break, end
 Gencount = Gencount + 1 ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end         %%%%%   End of 'While' Loop  %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


