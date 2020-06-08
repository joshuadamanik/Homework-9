%-----------------------ParamMaxX---------------------------------------------------
% Particle Swarm Optimization
%--------------------------------------------------------------------------
%
%% Version 1.0
% June 8, 2020
%
% Algorithm from https://en.wikipedia.org/wiki/Particle_swarm_optimization
% Modified from cealm_v20.m
%
% Compiled by:
% Joshua Julian Damanik (20194701)
% Aerospace Enginnering
% Korea Advanced Institute of Science and Technology (KAIST)

%---------- initial bounds of Lagrange multipliers

   NumParamY = NumIneq + NumEq;
   
   if(NumIneq ~= 0)
       for k=1:NumIneq
          ParamMaxY(k) = 0.01;
          ParamMinY(k) = 0;     
       end
   end    
   if(NumEq ~= 0)
       for k=1:NumEq
          ParamMaxY(NumIneq+k) = 0.01;
          ParamMinY(NumIneq+k) = -0.01;     
       end
   end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Gencount 	 = 1 ;  	                  % Generation counter

   InitFlag      = 0 ;

   IntervalX     = ParamMaxX - ParamMinX;
   IntervalY     = ParamMaxY - ParamMinY;
   
   VelocityMinX  = -IntervalX;
   VelocityMaxX  = IntervalX;
   
   VelocityMinY  = -IntervalY;
   VelocityMaxY  = IntervalY;

   OffspringPopX    = zeros(NumOffspringX,NumParamX);
   OffspringPopY    = zeros(NumOffspringY,NumParamY);
   
   VelocityPopX  = zeros(NumOffspringX,NumParamX);
   VelocityPopY  = zeros(NumOffspringY,NumParamY);
   
   BestPopX      = zeros(NumOffspringX,NumParamX);
   BestPopY      = zeros(NumOffspringY,NumParamY);

   Tol           = Tolerance*ones(1,NumParamY);              

   PrevF          = 1000000000000000.;
   PrevV          = 1000000000000000.;

%---------------------------------------------------
% Initial Population: 
%    Uniform distrubution within the search space
%---------------------------------------------------
  
%   randn('seed',sum(100*clock));                % Used 'seed' instead of 'state' for 
%   randn('seed',sum(100*clock));               % compatability of Matcom3 to Matlab4  

  for i=1:NumOffspringX  
     for j=1:NumParamX
         OffspringPopX(i,j) = rand*(ParamMaxX(j)-ParamMinX(j)) + ParamMinX(j);
         BestPopX(i,j) = OffspringPopX(i,j);
         VelocityPopX(i,j) = (rand*(VelocityMaxX(j)-VelocityMinX(j)) + VelocityMinX(j));
%          OffspringSigmaX(i,j)= VarIntervalX(j);
     end
  end

  for i=1:NumOffspringY   
     for j=1:NumParamY
        OffspringPopY(i,j) = rand*(ParamMaxY(j)-ParamMinY(j)) + ParamMinY(j);
        BestPopY(i,j) = OffspringPopY(i,j);
        VelocityPopY(i,j) = (rand*(VelocityMaxY(j)-VelocityMinY(j)) + VelocityMinY(j));
%         OffspringSigmaY(i,j)= VarIntervalY(j);
     end
  end   
  
SX = ceil(rand*NumOffspringX);
SY = ceil(rand*NumOffspringY);

BestParamX = BestPopX(SX,:);
BestParamY = BestPopY(SY,:);

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
    RP = rand();
    RG = rand();
    for j=1:NumParamX
        VelocityPopX(i,j) = omega * VelocityPopX(i,j)...
            + phiP * RP * (BestPopX(i,j) - OffspringPopX(i,j))...
            + phiG * RG * (BestParamX(j) - OffspringPopX(i,j));
%         VelocityPopX(i,j) = min(max(VelocityPopX(i,j),VelocityMinX(j)),VelocityMaxX(j));
        OffspringPopX(i,j) = OffspringPopX(i,j) + VelocityPopX(i,j);
    end
end

for i=1:NumOffspringY
    RP = rand();
    RG = rand();
    for j=1:NumParamY
        VelocityPopY(i,j) = omega * VelocityPopY(i,j)...
            + phiP * RP * (BestPopY(i,j) - OffspringPopY(i,j))...
            + phiG * RG * (BestParamY(j) - OffspringPopY(i,j));
%         VelocityPopY(i,j) = min(max(VelocityPopY(i,j),VelocityMinY(j)),VelocityMaxY(j));
        OffspringPopY(i,j) = OffspringPopY(i,j) + VelocityPopY(i,j);
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
         OffspringPopX(i,j) = OffspringPopX(j) - rand*min(violate1,IntervalX(j)) ;
      end
      if violate2 >0 
         OffspringPopX(i,j) = OffspringPopX(j) + rand*min(violate2,IntervalX(j)) ;
      end
   end
end

%-------------- lower bound check only for inequality constraints -------------- 

for i=1:NumOffspringY
   for j=1:NumIneq                        
      violate1 = OffspringPopY(i,j)-ParamMaxY(j);
      violate2 = ParamMinY(j)-OffspringPopY(i,j);
%       if violate1 >0 
%          OffspringPopY(i,j) = OffspringPopY(j) - rand*min(violate1,IntervalY(j)) ;
%       end
      if violate2 >0 
         OffspringPopY(i,j) = OffspringPopY(j)+rand*min(violate2,IntervalY(j)) ;
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
end

if exist('CostBestX','var') == 0
    BestPopX = OffspringPopX;
    BestPopY = OffspringPopY;
    
    BestCostPopX = maxcost;
    BestCostPopY = mincost;
    
    [BestCostX,bestindexx] = min(BestCostPopX);
    [BestCostY,bestindexy] = max(BestCostPopY);
    
    BestParamX = BestPopX(bestindexx,:);
    BestParamY = BestPopY(bestindexy,:);
end

for i=NumOffspringX
    if (maxcost(i) < BestCostPopX(i))
        BestPopX(i,:) = OffspringPopX(i,:);
        BestCostPopX(i) = maxcost(i);
    end
    if BestCostPopX(i) < BestCostX
        BestParamX = BestPopX(i,:);
        BestCostX = BestCostPopX(i);
    end
end

% [BestCostY, bestindexy] = max(CostOffspringY(:,bestindexx));
% BestParamY = OffspringPopY(bestindexy,:);

for i=NumOffspringY
    if (mincost(i) > BestCostPopY(i))
        BestPopY(i,:) = OffspringPopY(i,:);
        BestCostPopY(i) = mincost(i);
    end
    if BestCostPopY(i) > BestCostY
        BestParamY = BestPopY(i,:);
        BestCostY = BestCostPopY(i);
    end
end

%--------------------------------------------------
% Output to Display & File
%--------------------------------------------------


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
%  check_transient;
 graph_data;
%  if done == true, break, end
 Gencount = Gencount + 1 ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end         %%%%%   End of 'While' Loop  %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


