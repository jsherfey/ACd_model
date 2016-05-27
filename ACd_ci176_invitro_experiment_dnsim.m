% Model: ci176 (Type2 without h-current; "blue cells")
clear all % need this to clear persistent variable in holdstim()
% addpath(genpath('/home/jason/models/dnsim/Tallie-ACd-cells/dependencies'));
% addpath(genpath('~/code/dnsim'));
% cd(fileparts(mfilename('fullpath')));
cd /home/jason/models/dnsim/Tallie-ACd-cells/distributions/base
global BIOSIMROOT; BIOSIMROOT='~/code/dnsim';
% simulation parameters
dsfact=10; solver='euler'; usecoder=1; dt=.01;

%{

TARGET:
    AHPtime: 62.3333
    spikewidth: 1.68
    spikefreq: 19.6104
    threshfreq: 0.891424
    Vrest: -73.6875

MODEL:
AHPtime: 104
spikewidth: 2.46667
spikefreq: 19.7915
threshfreq: 0.902772
Vrest: -75.0811

%}

% INPUTS
if 0
  % tonic depolarization
  stepstim=0; depol=0; tspan=[0 30000]; % tspan=[0 30000];
  bltime=0; nsections=0; nsteps=0;
  noise=.15; tonicstim=.5;
  noise=.1; tonicstim=.505;
  fname=sprintf('ci176_tonic_noise%g',noise);
else
  % steps
  tonicstim=0; tspan=[0 60000]; depol=1;
  bltime=500; nsections=5; nsteps=5;
  stepstim=5; % 6
  noise=.06; % .15
  fname=sprintf('ci176_steps_noise%g',noise);
end

% additional step protocol parameters
lfac=1; stepreset=0;
isi=lfac*1000; steptime=lfac*400;
ramprate=1; membranearea=3000; tonictime=tspan(2); 
onset=50; % onset [ms] of tonic input
%I=getStepProtocolStim(.01,isi,nsteps,steptime,stepstim,membranearea,nsections,tonictime,bltime,tspan,ramprate,depol,stepreset); figure; plot(I)

% biophysical parameters
ECa=126.1; EK=-80; ENa=55; Cm=1;
  gnaf = 65;%50; % 50.25
  wbgNa=25;%25;
  gkdr = 4;%3.5; % 6
  gAHP = .01; % .01
  gnap = .0005;
  gks = 0;
  taurca = 100;%28.5714;
  gcat = .001;%.0025;
  gcan = 0; % .0056;
  gkca = .6;%.2; % 5
  taumin = 0;
  gh = 0; % .002
  eh = -10;
  gpas = 0.04; % .017
  epas = -75; % -66

mechanisms = {'cadyn','cat','kdr','AHP','iks','kca','naf','wbNa','stim','nap','pas','randn','iStepProtocol'};
%mechanisms = {'kdr','naf','pas','cadyn','cat','kca','randn','iStepProtocol','stim'};
%mechanisms = {'cadyn','cat','kca','kdr','AHP','naf','stim','nap','iks','h','pas','randn','iStepProtocol'};

spec=[];
spec.nodes(1).label = 'ci176';
spec.nodes(1).multiplicity = 1;
spec.nodes(1).dynamics = {'v''=current/c'};
spec.nodes(1).mechanisms = mechanisms;
spec.nodes(1).parameters = {'depol',depol,...
  'c',Cm,'v_IC',-65,'stim',tonicstim,'noise',noise,'onset',onset,'stepreset',stepreset,...
  'isi',isi,'nsteps',nsteps,'steptime',steptime,'stepsize',stepstim,'nsections',nsections,'tonictime',tonictime,'membranearea',membranearea,'bltime',bltime,'ramprate',ramprate,...
  'gnaf',gnaf,'gnap',gnap,'taurca',taurca,'gcat',gcat,'gcan',gcan,...
  'gkca',gkca,'taumin',taumin,'gkdr',gkdr,'gks',gks,'gAHP',gAHP,'wbgNa',wbgNa,...
  'gh',gh,'gpas',gpas,'epas',epas,'eh',eh,'ek',EK,'eca',ECa,'ena',ENa};

data = runsim(spec,'timelimits',tspan,'dt',dt,'dsfact',dsfact,'coder',usecoder,'SOLVER',solver);
AnalyzeStepProtocol; 
% plotv(data,spec,'varlabel','v');
% plotpow(data,spec,'NFFT',4096,'varlabel','v');
%plotpow(sim_data,spec,'NFFT',4096,'varlabel','v','FreqRange',[1 100]);

return

outpath='/home/jason/models/dynasim/ACC_simulations/paper1_results';
file=fullfile(outpath,[fname '_winner']);
print(gcf,[file '.jpg'],'-djpeg'); print(gcf,[file '.eps'],'-depsc'); saveas(gcf,[file '.fig']); plot2svg([file '.svg'],gcf);
save([file '.mat']);

% re-run without steps
if 0
  args={...
  'c',Cm,'v_IC',-65,'stim',tonicstim,'noise',noise,'onset',onset,'stepreset',stepreset,...
    'isi',isi,'nsteps',nsteps,'steptime',steptime,'stepsize',stepstim,'nsections',nsections,'tonictime',tonictime,'membranearea',membranearea,'bltime',bltime,'ramprate',ramprate,...
    'gnaf',gnaf,'gnap',gnap,'taurca',taurca,'gcat',gcat,'gcan',gcan,...
    'gkca',gkca,'taumin',taumin,'gkdr',gkdr,'gks',gks,'gAHP',gAHP,...
    'gh',gh,'gpas',gpas,'epas',epas,'eh',eh,'ek',EK,'eca',ECa,'ena',ENa};
  tonicstim=.9; tspan=[0 1000];
  spec.nodes(1).mechanisms = {'cadyn','cat','kca','kdr','AHP','naf','stim','nap','iks','h','pas','randn'};
  spec.nodes(1).parameters = {args{:},'depol',0,'stim',tonicstim};
  data = runsim(spec,'timelimits',tspan,'dt',dt,'dsfact',dsfact,'coder',usecoder,'SOLVER',solver);
  %AnalyzeStepProtocol; %plotv(data,spec,'varlabel','v');
  plotv(data,spec,'varlabel','v');
end

% Type2:          {'cadyn','cat','kca','kdr','naf','stim','nap','iks','h','pas','noise','iStepProtocol'};
% PFC_models:     {'naf','nap','cadyn','cat','can','kca','kdr','iks','AHP','h','pas','stim','randn','iStepProtocol'};
% PFC_models_min: {'naf','E_Nap','wbNa','cadyn','can','kca','kdr','M','AHP','h','pas','stim','randn','iStepProtocol'};

%{
TARGET:
    AHPtime: 62.3333
    spikewidth: 1.68
    spikefreq: 19.6104
    threshfreq: 0.891424
    Vrest: -73.6875
  [Idep;fdep]
    100.0000  200.0000  300.0000  400.0000  500.0000
           0         0    8.8652   11.4155   15.2439
  [Idep/100;fdep0]
    100.0000  200.0000  300.0000  400.0000  500.0000
           0         0   52.0833   94.3396   34.2466
  [Idep/100;fdepSS]
    100.0000  200.0000  300.0000  400.0000  500.0000
           0         0    7.3421    9.9602   12.0482
  [I/100;v]
   -500.0000 -400.0000 -300.0000 -200.0000 -100.0000  100.0000
    -95.1248  -92.5744  -88.2141  -83.7677  -78.7946  -67.6256
%}