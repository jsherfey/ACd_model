clear all
cd /home/jason/models/dynasim/ACC_simulations/paper1_results
Ne=80; Ni=.25*Ne;
f1=[0:7.5:60]; f2=[0:7.5:60]; acAMPA1e=4500; enoise=acAMPA1e; 
tauGABA=5;
gAMPAee=[.05]; % E->E
gNMDAee=[.2];
gAMPAei=1;
tspan=[0 5000]; dt=.01; solver='rk1'; compile_flag=1; downsample_factor=10;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coupled_type={[0 0],[0 1],[1 1]}; % gEE:(between,within)
G=[2 3]; % coupled types to process
hetdegrees=[0 .05 .1 .2 1];       % fraction of data-constrained covariance
I=[1 2 3 4 5]; % hetdegrees to process; eg) het(1)=0 vs het(2)=.05
I=[1 4];
colors='brcgk'; % hom/het colors
reps=1:5; % repetitions to process (max: reps=1:5)
h1=figure('position',[145 100 1520 820]); nrows=2; ncols=length(coupled_type);
h2=figure('position',[145 100 1520 820]);
plottype=3;
  % 1. scatter plot
  % 2. means per df
  % 3. mean+/-std per df
dN=[];
for g=G%2:length(coupled_type)
  betweenblockee=coupled_type{g}(1);
  withinblockee=coupled_type{g}(2);
  for i=I %1:length(hetdegrees)
    hetdegree=hetdegrees(i);
    for j=1:length(reps)
      rep=reps(j);
      study_dir_root=sprintf('sweeps_het-vs-hom_ac%g_EEwithin%g-%g_EEbetween%g-%g_gEEnmda%g-%g_gEEampa%g-%g_gEI%g_%g-%gms',acAMPA1e,withinblockee(1),withinblockee(end),betweenblockee(1),betweenblockee(end),gNMDAee(1),gNMDAee(end),gAMPAee(1),gAMPAee(end),gAMPAei,tspan);
      if hetdegree==0
        study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HOM_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,rep);
      else
        study_dir=sprintf('%sNe%g_f_%g-%gHz-vs-%g-%gHz_tauI%gms_HET%g_%g',study_dir_root,Ne,f1(1),f1(end),f2(1),f2(end),tauGABA,hetdegree,rep);
      end
      color=colors(i);
      statfile=[study_dir '_stats.mat'];
      if ~exist(statfile,'file'), continue; end
      load(statfile,'stats');
      % calculate metrics (competition,coherence)
      F1=[stats.E_fAMPA1]; uf1=unique(F1); nf1=length(uf1);
      F2=[stats.E_fAMPA2]; uf2=unique(F2); nf2=length(uf2);
      f0i=(F1==0)|(F2==0); % comparison to nonrhythmic or not
      n1=arrayfun(@(x)x.pairs.Nspikes1,stats);
      n2=arrayfun(@(x)x.pairs.Nspikes2,stats);
      dn=n1-n2;
      dn=(n1-n2)./max(n1,n2);
      df=F1-F2;
      c12=arrayfun(@(x)x.pairs.xcsum_pops,stats);
      if isempty(dN)
        udf=unique(F1-F2);
        ndf=length(udf);
        dN =cell(ndf,length(coupled_type),length(hetdegrees)); dN0=dN;
        C12=cell(ndf,length(coupled_type),length(hetdegrees)); C120=C12;
      end
      for k=1:ndf
        ix=(df==udf(k))&(~f0i);
        dN{k,g,i}=[dN{k,g,i} dn(ix)']; % mean(dn(ix));
        C12{k,g,i}=[C12{k,g,i} c12(ix)'];
        ix=(df==udf(k))&(f0i);
        dN0{k,g,i}=[dN0{k,g,i} dn(ix)']; % mean(dn(ix));
        C120{k,g,i}=[C120{k,g,i} c12(ix)'];
      end
      % plot metrics
      if plottype==1 % SCATTER PLOTS
        subplot(nrows,ncols,g); % competition
        plot(df(~f0i),abs(dn(~f0i)),[color 'o'],df(f0i),abs(dn(f0i)),[color 'x']); hold on
        xlabel('df'); ylabel('|dn|'); line(xlim,[0 0]); line([0 0],ylim);
        subplot(nrows,ncols,g+ncols); % coherence 
        plot(df(~f0i),c12(~f0i),[color 'o'],df(f0i),c12(f0i),[color 'x']); hold on
        xlabel('df'); ylabel('c12'); ylim([0 1]); line([0 0],ylim); 
      end
    end
    if plottype==2 % MEANS per df
      T=1;%(tspan(2)-tspan(1))/1000/Ne; % convert # spikes in pop to <rates>
      dnabsmu=cellfun(@(x)mean(abs(x/T)),dN(:,g,i));
      dn0absmu=cellfun(@(x)mean(abs(x/T)),dN0(:,g,i));
      c12mu=cellfun(@mean,C12(:,g,i));
      c120mu=cellfun(@mean,C120(:,g,i));
      subplot(nrows,ncols,g); % competition
      plot(udf,dnabsmu,[color 'o-'],udf,dn0absmu,[color 'x-']); ylim([0 1]); hold on
      xlabel('df'); ylabel('<%dn> [Hz]'); line(xlim,[0 0],'color','k'); line([0 0],ylim,'color','k');
      subplot(nrows,ncols,g+ncols); % coherence 
      plot(udf,c12mu,[color 'o-'],udf,c120mu,[color 'x-']); ylim([0 1]); hold on
      xlabel('df'); ylabel('<c12>'); ylim([0 1]); line([0 0],ylim,'color','k'); 
    elseif plottype==3 % MEAN+/-STD per df, + markers for significance (ttest2)
      % mu
      dnabsmu=cellfun(@(x)mean(abs(x)),dN(:,g,i));    dnabsmu(isnan(dnabsmu))=0;
      dn0absmu=cellfun(@(x)mean(abs(x)),dN0(:,g,i));  dn0absmu(isnan(dn0absmu))=0;
      c12mu=cellfun(@mean,C12(:,g,i));                c12mu(isnan(c12mu))=0;
      c120mu=cellfun(@mean,C120(:,g,i));              c120mu(isnan(c120mu))=0;
      % sd
      dnabssd=cellfun(@(x)std(abs(x)),dN(:,g,i));     dnabssd(isnan(dnabssd))=0;
      dn0abssd=cellfun(@(x)std(abs(x)),dN0(:,g,i));   dn0abssd(isnan(dn0abssd))=0;
      c12sd=cellfun(@std,C12(:,g,i));                 c12sd(isnan(c12sd))=0;
      c120sd=cellfun(@std,C120(:,g,i));               c120sd(isnan(c120sd))=0;
      % se
      dnabsse=dnabssd./sqrt(max(1,cellfun(@numel,dN(:,g,i))));
      dn0absse=dn0abssd./sqrt(max(1,cellfun(@numel,dN0(:,g,i))));
      c12se=c12sd./sqrt(max(1,cellfun(@numel,C12(:,g,i))));
      c120se=c120sd./sqrt(max(1,cellfun(@numel,C120(:,g,i))));
      % plots
      figure(h1);
      subplot(nrows,ncols,g); % competition
      sel=1:length(udf);
      sel(udf==0)=[];
      sel=sel(2:end-1);
      h=plot_CI(udf(sel),dnabsmu(sel)*100,dnabsse(sel)*100,color,'o'); hold on; ylim([0 .5]*100);
      subplot(nrows,ncols,g+ncols); % coherence 
      h=plot_CI(udf(sel),c12mu(sel),c12se(sel),color,'o'); hold on; ylim([0 .4]);
      figure(h2);
      subplot(nrows,ncols,g); % competition
      sel=1:length(udf);
      sel(udf==0)=[];
      h=plot_CI(udf(sel),dn0absmu(sel)*100,dn0absse(sel)*100,color,'o'); hold on; ylim([0 1]*100);
      subplot(nrows,ncols,g+ncols); % coherence 
      h=plot_CI(udf(sel),c120mu(sel),c120se(sel),color,'o'); hold on; ylim([0 .4]);
    end
    
  end
end

%% significance
g=2; I=[1 4]; % compare het=0 vs het=20% for coupled(0,1)
% competition: f1 vs f2
X=dN(:,g,I(1));
Y=dN(:,g,I(2));
P=nan(1,ndf);
for k=1:ndf
  [h,p]=ttest2(abs(X{k}),abs(Y{k}));
  if ~isempty(p)
    P(k)=p;
  end
end
sigf=udf(P<.05);
figure(h1); subplot(nrows,ncols,g); plot(sigf,45,'r*');

% competition: f1 vs 0
X=dN0(:,g,I(1));
Y=dN0(:,g,I(2));
P=nan(1,ndf);
for k=1:ndf
  [h,p]=ttest2(abs(X{k}),abs(Y{k}));
  if ~isempty(p)
    P(k)=p;
  end
end
sigf=udf(P<.05);
figure(h2); subplot(nrows,ncols,g); plot(sigf,45,'r*');

% coherence: f1 vs f2
X=C12(:,g,I(1));
Y=C12(:,g,I(2));
P=nan(1,ndf);
for k=1:ndf
  [h,p]=ttest2(abs(X{k}),abs(Y{k}));
  if ~isempty(p)
    P(k)=p;
  end
end
sigf=udf(P<.05);
figure(h1); subplot(nrows,ncols,g+ncols); plot(sigf,.35,'r*');

% coherence: f1 vs 0
X=C120(:,g,I(1));
Y=C120(:,g,I(2));
P=nan(1,ndf);
for k=1:ndf
  [h,p]=ttest2(abs(X{k}),abs(Y{k}));
  if ~isempty(p)
    P(k)=p;
  end
end
sigf=udf(P<.05);
figure(h2); subplot(nrows,ncols,g+ncols); plot(sigf,.35,'r*');

return

figure(h1);
filename='coupled-decoupled-vs-DF_hom-vs-het_without-0Hz_with-significance.svg';
plot2svg(filename,gcf);

figure(h2);
filename='coupled-decoupled-vs-DF_hom-vs-het_with-0Hz_with-significance.svg';
plot2svg(filename,gcf);


