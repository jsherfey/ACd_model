
addpath /projectnb/crc-nak/sherfey/projects/ACC_simulations/ACC_models % add path to getEE()
cd /projectnb/crc-nak/sherfey/projects/ACC_simulations/ACC_models
ipfile='rat-ACd-IPs/rat-ACd-model_cell-intrinsic-properties.mat';

%% Model specification

Ne=80; % 10,20,40,80,100
Ni=.25*Ne;

% E-cell mechanisms
mechanisms={'AHP','cadyn','can','cat','h','iks','kdr','kca','naf','nap','pas'};

% load model parameters fit to single cell intrinsic property data
load(ipfile,'simdata','expdata');

% set baseline E-cell parameter values (median) for homogeneous network
parameters={};
for i=1:length(simdata.parameters)
  param=simdata.parameters{i};
  value=median(simdata.(param));
  parameters{end+1}=param;
  parameters{end+1}=value*ones(1,Ne); % homogeneous across population
end
% account for default-value discrepancies b/w DNSim [PP13] and DynaSim [DS02] implementations
aAHP_scale=1; % 1, 1e6
CAF=600*4.6452; % 600, 600*4.6452
cainf = 5e-5; % .05, 5e-5
parameters{end+1}='aAHP_scale';
parameters{end+1}=aAHP_scale;
parameters{end+1}='CAF';
parameters{end+1}=CAF;
parameters{end+1}='cainf';
parameters{end+1}=cainf;
parameters{end+1}='IC_noise';
% set initial condition noise level (IC_noise*rand offsets each cell's gating variable IC)
IC_noise=.25;
parameters{end+1}=IC_noise;

% collect values for parameters to make heterogeneous in E-cell population
het={'gAHP','gcan','gcat','gh','gkca','gks','gkdr','gnaf','gpas'};
params=[];
for i=1:numel(het)
  params=cat(2,params,simdata.(het{i}));
end
% parameters of the multivariate normal distribution
MU = median(params,1);
SIGMA = cov(params);

% define random seeds for all simulations (use same for hom and het)
num_seeds=1e3; % <= max_num_realizations
seedfile=sprintf('%gseeds.mat',num_seeds);
if exist(seedfile,'file')
  load(seedfile,'seeds');
else
  seeds=zeros(num_seeds,1);
  for i=1:num_seeds
    rng('shuffle');
    seeds(i)=getfield(rng,'Seed');
  end
  % save seeds to file (for reproducibility)
  save(seedfile,'seeds');
end

% generic inputs and state equations
assembly_size=Ne/2;
Ke1=zeros(1,Ne); Ke1(1:assembly_size)=1;    % input kernel for first E-cell assembly
Ke2=zeros(1,Ne); Ke2(assembly_size+1:Ne)=1; % input kernel for all other E-cell assemblies
E_input_def={'input(V)=-gAMPA.*(s1(k,:)+s2(k,:)).*(X-EAMPA); monitor input; EAMPA=0; gAMPA=0; onset=50; offset=inf;';
     sprintf('s1=getPoissonGating(baseline/2,dcAMPA1,acAMPA1,fAMPA1,phiAMPA1,onset,offset,tauAMPA,T,Npop,%s,kick,ramp_dc_flag1,ramp_ac_flag1);',toString(Ke1));
     sprintf('s2=getPoissonGating(baseline/2,dcAMPA2,acAMPA2,fAMPA2,phiAMPA2,onset,offset,tauAMPA,T,Npop,%s,kick,ramp_dc_flag2,ramp_ac_flag2);',toString(Ke2));
             'baseline=0; dcAMPA1=0; acAMPA1=0; fAMPA1=0; phiAMPA1=0; tauAMPA=2; kick=1; ramp_dc_flag1=0; ramp_ac_flag1=0; dcAMPA2=0; acAMPA2=0; fAMPA2=0; phiAMPA2=0; ramp_dc_flag2=0; ramp_ac_flag2=0;';
            };
I_input_def={'input(V)=-gAMPA.*s1(k,:).*(X-EAMPA); monitor input; EAMPA=0; gAMPA=0; onset=50; offset=inf;';
             's1=getPoissonGating(baseline,dcAMPA1,acAMPA1,fAMPA1,phiAMPA1,onset,offset,tauAMPA,T,Npop,ones(1,Npop),kick,ramp_dc_flag1,ramp_ac_flag1);';
             'baseline=0; dcAMPA1=0; acAMPA1=0; fAMPA1=0; phiAMPA1=0; tauAMPA=2; kick=1; ramp_dc_flag1=0; ramp_ac_flag1=0;';
            };
state_equations='dV/dt=(@current+input(V)+Iapp*(t>onset&t<offset))./Cm; Cm=1; Iapp=0; V(0)=-65;';

% input parameters
onset=100;   % ms, stimulus start time
offset=inf;  % ms, stimulus stop time
gAMPAe=.001; % mS/cm2
enoise=4500; % Hz
dcAMPA1e=0;     dcAMPA2e=dcAMPA1e; % Hz
acAMPA1e=4500;  acAMPA2e=acAMPA1e; % Hz
fAMPA1e=0;      fAMPA2e=0;
E_input_parameters={'dcAMPA1',dcAMPA1e,'dcAMPA2',dcAMPA2e,'acAMPA1',acAMPA1e,'acAMPA2',acAMPA2e,'fAMPA1',fAMPA1e,'fAMPA2',fAMPA2e,'gAMPA',gAMPAe,'baseline',enoise,'tauAMPA',tauAMPA,'onset',onset,'offset',offset};
gAMPAi=0;
inoise=0;
dcAMPA1i=0;
acAMPA1i=0;
fAMPA1i=0;
I_input_parameters={'dcAMPA1',dcAMPA1i,'acAMPA1',acAMPA1i,'fAMPA1',fAMPA1i,'gAMPA',gAMPAi,'baseline',inoise,'tauAMPA',tauAMPA,'onset',onset,'offset',offset};

% connectivity parameters and kernels
withinblockee=1;  % scaling for connectivity kernels within each assembly
betweenblockee=0; % scaling for connectivity kernels between assemblies
gAMPAee=.05; % E->E (<=.1)
gNMDAee=.2; % E->E (.1-.2)
gAMPAie=.1; % I->E (.1)
gAMPAii=1; % I->I
gAMPAei=1; % E->I (.1)

% define and normalize kernels by number of presynaptic connections
Kii=ones(Ni)-eye(Ni);
Kei=ones(Ne,Ni);
Kie=ones(Ni,Ne);
Kei=Kei./repmat(max(1,sum(Kei,1)),[size(Kei,1) 1]);
Kie=Kie./repmat(max(1,sum(Kie,1)),[size(Kie,1) 1]);
Kii=Kii./repmat(max(1,sum(Kii,1)),[size(Kii,1) 1]);
% Kee: E->E kernels are set by getEE() called in iAMPAee.mech and iNMDAee.mech

% synaptic time constants
tauGABA=5; % ms, inhibition decay time constant
tauAMPA=2; 
tauNMDA=95;

% NETWORK SPECIFICATION
spec=[];
% Populations:
% E-cells (heterogeneous parameters)
spec.populations(1).name='E';
spec.populations(1).size=Ne;
spec.populations(1).equations=[state_equations E_input_def{:}];
spec.populations(1).mechanism_list=mechanisms;
spec.populations(1).parameters={parameters{:},E_input_parameters{:}};
% Wang-Buzsaki interneuron model, 1996
spec.populations(2).name='I';
spec.populations(2).size=Ni;
spec.populations(2).equations=[state_equations I_input_def{:}];
spec.populations(2).mechanism_list={'WB96FSiNa','WB96FSiK','WB96FSileak'};
spec.populations(2).parameters={'Cm',1,'Eleak',-65,'gleak',.1,'gNa',35,'gK',9,I_input_parameters{:}};    
% Connections:
% E->E (AMPA,NMDA): 50% assemblies
spec.connections(1).direction='E->E';
spec.connections(1).mechanism_list={'iAMPAee','iNMDAee'};
spec.connections(1).parameters={'gAMPA',gAMPAee,'gNMDA',gNMDAee,'tauAMPA',tauAMPA,'tauNMDA',tauNMDA,'gw',withinblockee,'gb',betweenblockee};
% E->I (AMPA,NMDA): all-to-all
spec.connections(2).direction='E->I';
spec.connections(2).mechanism_list={'iAMPA'};
spec.connections(2).parameters={'gAMPA',gAMPAei,'netcon',Kei,'tauAMPA',tauAMPA};
% I->E (GABA): all-to-all
spec.connections(3).direction='I->E';
spec.connections(3).mechanism_list={'iGABA'};
spec.connections(3).parameters={'gGABA',gAMPAie,'netcon',Kie,'tauGABA',tauGABA};
% I->I (GABA): all-to-all
spec.connections(4).direction='I->I';
spec.connections(4).mechanism_list={'iGABA'};
spec.connections(4).parameters={'gGABA',gAMPAii,'netcon',Kii,'tauGABA',tauGABA};

% store baseline model
base=spec;

%% Simulations

% parameters to vary across simulations
f1=(0:7.5:60); 
f2=(0:7.5:60);
betweenblockee=0; % 0,1 kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
withinblockee=1; % 0,1
vary={'E->E','gw',withinblockee;'E->E','gb',betweenblockee;'E','fAMPA1',f1;'E','fAMPA2',f2;'E','(acAMPA1,acAMPA2)',acAMPA1e;'E','(dcAMPA1,dcAMPA2)',dcAMPA1e;'(I->E,I->I)','tauGABA',tauGABA};

% simulator options
tspan=[0 5000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
save_data_flag=0; cluster_flag=1; memory_limit='32G'; sims_per_job=10;
analysis_functions=@CalcSpikeSync;
analysis_options={'ROI_pairs',{'E_V',[0 .5],'E_V',[.5 1]}};
plot_functions={@PlotData,@PlotData};
plot_options={{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
simulator_options={'save_data_flag',save_data_flag,'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1,'cluster_flag',cluster_flag,'plot_functions',plot_functions,'plot_options',plot_options,'memory_limit',memory_limit,'downsample_factor',downsample_factor,'analysis_functions',analysis_functions,'analysis_options',analysis_options,'sims_per_job',sims_per_job};

% limit the number of cores to <= scc limit
maxNumCompThreads(4);

% run simulations
for rep=1:5 % loop over repetitions (each with a unique random seed)
  
  % reset model
  spec=base; 
  
  % HOMOGENEOUS ASSEMBLIES
  study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
  hom_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
  [data,hom_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',hom_study_dir,'vary',vary,simulator_options{:});

  % HETEROGENEOUS ASSEMBLIES
  hetdegree=.2;
  rng(seeds(rep)); % set the random seed for drawing parameter values
  values = mvnrnd(MU,hetdegree*SIGMA,Ne/2); % parameter values for one realization of E-cell population
  values = cat(1,values,values);  % copy for each assembly
  values = max(0,values);         % turn off currents w/ negative max conductance
  % update model specification with heterogeneous parameter values
  mods={};
  for i=1:numel(het)
    mods=cat(1,mods,{'E',het{i},values(:,i)'});
  end
  spec=ApplyModifications(spec,mods);
  study_dir_root=sprintf('studies/sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
  het_study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
  [data,het_studyinfo]=SimulateModel(spec,'random_seed',seeds(rep),'study_dir',het_study_dir,'vary',vary,simulator_options{:});

end

% unix(sprintf('cat %s/pbsout/sim_job1.out',hom_studyinfo.simulations(1).batch_dir));
% unix(sprintf('cat %s/pbsout/sim_job1.out',het_studyinfo.simulations(1).batch_dir));
