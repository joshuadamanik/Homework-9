%
%
%%%  Main Program of CEALM & PSO

  clear; clc; close all;

%------ Problem to be solved -------------------
%  fprintf('\n prob_m10 \n' ) ;
  prob_m01;                                 % defines the problem to be solved1
 
%--------- Monte-Carlo computation -----------------------------------------------
   nrun          = 10 ;                   % no. of MC computation
   MaxGen        = 6000;               % Maximum Generation 

   nPrint = MaxGen * 10;     % print every nPrint generations

%------ Output control -----------------------------------------
  outfile = fopen ('out_v20.dat','w');


%--------- Strategy -----------------------------------------------------

   method = 'pso';
%  pso   = Particle swarm optimization
%  cealm = Co-Evolution Augmented Lagrangian Method

   istrategy = 0;    
%  ONLY FOR CEALM. PSO ALWAYS USE SECURITY STRATEGY
%  0 = security strategy for both players (recommended)   
%  1 = man-to-man strategy for the parameter, which is the follower
%



%--------- Populations -----------------------------------------------------

   NumOffspringX  = 40;                 % CEALM Offspring & PSO Population
   NumParentX     = 8;                  % CEALM Parent Population
    
   NumOffspringY  = 40;                 % CEALM Offspring & PSO Population
   NumParentY     = 8;                  % CEALM Parent Population
 
%----------------------------------------------------------------------

  done = false;
  for irun=1:nrun
    fprintf(outfile,'\n Run=%3d ', irun) ;
    fprintf(outfile,'\n') ;
    if strcmp(method, 'pso')
        pso_v1
    else
        cealm_v20
    end
    
  end
  done = true;
  graph_data;

  fclose(outfile);

%--------- end -----------

