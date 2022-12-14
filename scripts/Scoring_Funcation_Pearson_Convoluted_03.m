close all;
ReadData = 1;           % Slow part of the script is reading the theoretical data. Set to 1 to read in data. 
if ReadData ==  1; 
    clear all;  ReadData = 1; 
end

%% ==================================================================================
% Order the simulations so that the grid-plot makes sense. 
 stepnr = [531 531 502 532 503 533 504 534 505 535 ...
                506 536 507 537 508 538 508 539 510 540 ...
                511 541 512 542 513 543 514 544 515 545 ...
                516 546 517 547 518 548 519 549 520 550 ...
                521 551 522 552 523 553 524 554 525 555 ...
                526 556 527 557 528 558 529 559 530 560];    
% This lists the simulations in the correct order for plotting.                                                                      
FitStart = 20; FitEnd = 350;              %FitStart = 20 in paper at this point.                   % Low q and high q limits used to fit experimental data. 
corrStart = 42; corrEnd = 348; 
ScaleCorr1 = 1.005;                 % best correlation between data but tends to stretch peak.   
q = [0:0.01:2]*17.7/18;         % Correct X-ray energy for the weighted average and not peak. 
KeepData = 100; 
BestScore = 0;                      % Comment out except for first time.
BestYet = 1; 
%% ===================================================================================
% Read in experimental data from Richard's analysis. 
        TEMP1 = load('/path/bRbasisRN.mat'); 
        Expt1 = TEMP1.Data2save; 
         Expt1(2:3,:) = Expt1(2:3,:);           
         % Expt1(3,:)  = Expt1(2,:);  % Comment to Fit second component; uncomment to fit first component. 
        q1 = ScaleCorr1*Expt1(1,:); 
% Read in Solvent Exclusion term from Solvent_Correction_Gromacs.m
        TEMP10 = load('/path/SolventExclusionFactor.mat'); 
        SolventExclusionFactor = TEMP10.SolventExclusionFactor2;   %  TempExpt = transpose(TempExpt);      
% Read in membrane correction from the BOG titration.         
        TEMP101 = load('/path/MembraneComponent2.mat'); 
        FlucMicel = 10^-8*TEMP101.Data2save; 
% Read in membrane correction from the micelle fluctuations.      
        TEMP102 = load('/path/MicelleFluctuationsFinalQuarter.mat'); 
        SVDmicelle = TEMP102.Data2save;     
%   FlucMicel = SVDmicelle*max(FlucMicel(2,:))/max(SVDmicelle(2,:));    
% Read in Undulator spectrum     
       TEMP103 = load('/path/U17.mat');    
       U17 = TEMP103.Data2save;   
       
%% ===================================================================================
% Read in results from Daniel's Gromacs analysis 
if ReadData == 1; 
for kk =  1:size(stepnr,2)    
        filename1 = ["/path/bR-grid5-noSCposres/input_fit/bR-190BOG-501/input_fit/bR-190BOG-rest-mem-rest-prot/bR-190BOG-rest-mem-rest-prot"] ; 
        filename2 = ['/path/bR-grid5-noSCposres/input_fit/bR-190BOG-' num2str(stepnr(kk)) '/input_fit/bR-190BOG-' num2str(stepnr(kk)) '-exc-mem-exc-prot/bR-190BOG-' num2str(stepnr(kk)) '-exc-mem-exc-prot'] ;
        filename3 = ['/path/bR-grid5-noSCposres/input_fit/bR-190BOG-' num2str(stepnr(kk)) '/input_fit/bR-190BOG-' num2str(stepnr(kk)) '-exc-mem-rest-prot/bR-190BOG-' num2str(stepnr(kk)) '-exc-mem-rest-prot'] ;  
        filename4 = ['/path/bR-grid5-noSCposres/input_fit/bR-190BOG-' num2str(stepnr(kk)) '/input_fit/bR-190BOG-' num2str(stepnr(kk)) '-rest-mem-exc-prot/bR-190BOG-' num2str(stepnr(kk)) '-rest-mem-exc-prot'] ;  
        filename = [filename1 filename2 filename3 filename4];
% Protein moving, membrane moving;   
        Load = [filename(1) filename(2)];
        DataDI = LoadDiffData(Load)*10^-8;    % Function LoadDiffData defined at bottom.
% Convolute these predictions with U17 undulator spectrum        
        for i = 1:size(DataDI,3)
            DataDI(6,:,i) = ConvU17(DataDI(6,:,i),q,U17);
        end              
%  Pick the conformations that we want to work with from the vacuum result.
        for i = 1:size(DataDI,3)
            test2(i,:) = [corr2(interp1(q,DataDI(6,:,i).*SolventExclusionFactor,q1(corrStart:corrEnd)),Expt1(3,corrStart:corrEnd)),i]; 
        end 
            [Xsorted,Ysorted] = sort(test2,1,'descend');
% Use this information to decide what sub-trajectories to keep
        for i=1:KeepData
            DSgood(i,:) = DataDI(6,:,Ysorted(i,1));
        end   
           dProtVac(kk,:) = mean(DSgood);      
% Load the cross-term between micelle and protein. 
        clear('DataDI')                                                               
        Load = [filename(1) filename(4)];
        DataDI = LoadDiffData(Load)*10^-8;   
% Convolute these predictions with U17 undulator spectrum        
        for i = 1:size(DataDI,3)
            DataDI(2,:,i) = ConvU17(DataDI(2,:,i),q,U17);
        end      
% Keep the information chosen above.         
        for i=1:KeepData
            DSgood(i,:) = DataDI(2,:,Ysorted(i,1));
        end   
            dPmMfgvac(kk,:) = mean(DSgood);          
% Load the complementary cross-term between micelle and protein.      
        clear('DataDI')
        Load = [filename(3) filename(2)];
        DataDI = LoadDiffData(Load)*10^-8; 
% Convolute these predictions with U17 undulator spectrum        
        for i = 1:size(DataDI,3)
            DataDI(2,:,i) = ConvU17(DataDI(2,:,i),q,U17);
        end      
% Keep the information chosen above.      
        for i=1:KeepData
            DSgood(i,:) = DataDI(2,:,Ysorted(i,1));
        end    
            dPmMfevac(kk,:) = mean(DSgood);               
end 
end
%% Evaluate the correlations on the vacuum term. 
        for i = 1:size(stepnr,2)   
          InterpTheory(i,:) = interp1(q,dProtVac(i,:).*SolventExclusionFactor,q1(corrStart:corrEnd)); 
          Correlations(i,:) = corr2(InterpTheory(i,:),Expt1(3,corrStart:corrEnd)); 
        end 
max(Correlations(1:size(stepnr,2),:));
[X1 X2] = max(Correlations); 

%% ===============================================================================================================
% The following loop does the fitting of the data & searchers for optimal
% parameters 
for ll = 1:21
          Bfactor = 50; % Best for State 1 & 2
          % Bfactor = 20 + (ll-11);
           DampWidth = 0.129;   % Best for State  1 & 2 
           % DampWidth = 0.129 + (ll-11)/1000;     % Allows micell damping 
            % BOGcorr = -173; % Best for State 1
            BOGcorr = -100; % Best for State 2
            % BOGcorr = -100 + (ll-11);    % Correction for deviations from 190 BOG simulation. 
for kk = 1:size(stepnr,2)
%% ==========================================================================
 % Prepare the matrices we use to score the fit.
   % LowDamp = normpdf(q,0,DampWidth); LowDamp = LowDamp/max(LowDamp);  
     LowDamp = (1-erf((q-DampWidth)/(DampWidth/4)))/2;
     % LowDamp = q*0; % Done to get plot without this correction. 
     AllDamp = normpdf(q,0,Bfactor); AllDamp = AllDamp/max(AllDamp);  
     dProtMemVac(kk,:) = (dPmMfevac(kk,:) + dPmMfgvac(kk,:))/2;               % This is what is used to fit our data. 
     dMemVac(kk,:) = (dProtMemVac(kk,:) + BOGcorr*FlucMicel(2,:) - dProtVac(kk,:)).*LowDamp;                        
     dProtMemVac(kk,:) =  AllDamp.*(dProtVac(kk,:) + dMemVac(kk,:));      % Damp the micell correction.
%  Interpolate theoretical curves onto the data to allow fitting. 
    dProtMemVacInterp(:,kk) = interp1(q,dProtMemVac(kk,:),q1(FitStart:FitEnd)); 
    Contrast = transpose(interp1(q,SolventExclusionFactor,q1(FitStart:FitEnd)));
% Shorten the experimental vectors to the chosen fitting lengths. 
    ExptData2 = transpose(Expt1(2,FitStart:FitEnd)); 
    ExptData3 = transpose(Expt1(3,FitStart:FitEnd)); 
% Shorten the theoretical vectors accordingly    
%     dProtMemVacInterp2(:,kk) = dProtMemVacInterp(FitStart:FitEnd,kk);
% Minimize against experimental data without low-q damping of membrane.  
    weight = q1(FitStart:FitEnd).^2.*gaussmf(q1(FitStart:FitEnd),[0.4 -0]);  
    weight = transpose(weight/max(weight));
    % close all; plot(q1(FitStart:FitEnd),weight)
    % Make the fit optimizing the correlation function. 
                score2(kk) =  corr2(dProtMemVacInterp(:,kk).*Contrast,ExptData3);     
% Find the best scaling values for the most correlated componet.             
                    fun = @(x)(sum((x(1)*dProtMemVacInterp(:,kk).*Contrast.*weight - ExptData3.*weight).^2));               
                    x0 = [20 0.0005];  x1 = fminsearch(fun,x0);     
% Write the correctly scaled function for this iteraction for plotting. 
            score300 = sum((ExptData3.*weight).^2); 
            scoreComp1(kk) = sqrt(sum((x1(1)*dProtMemVacInterp(:,kk).*Contrast.*weight - ExptData3.*weight).^2))/sqrt(score300);
% Check to see if it is an improvement, and if so then save values. 
            if  scoreComp1(kk) < BestYet; 
                    BestYet = scoreComp1(kk); 
                    Best3(kk,:) = x1(1)*dProtMemVacInterp(:,kk).*Contrast;    
                    BestIndex = kk;       
                    BestScale = x1(1);
                    BestDampWidth = DampWidth; 
                    BestBfactor = Bfactor; 
                    BestBOGcorrection = BOGcorr;
            end  
end 
end
% Print to screen the optimal values. 
           BestDampWidth 
           % BestBfactor 
           BestBOGcorrection
           BestScale = 30; 
 % Summarise the results of scaling components 1 & 2 together.            
         ScoreTogether = [ 20 35.5 29.6; 
                                       25 35.9 22.4;
                                       27.5 33.9 22.6; 
                                       30 33.6 20.4; 
                                       32.5 35.6 20.1; 
                                       35 38.0 19.9;  
                                       40 42.0 19.3; 
                                       45 46.3 22.8; 
                                       50 51.6 17.9]; 
                                   
              close all; plot(ScoreTogether(:,1),sqrt(ScoreTogether(:,2).^2+ScoreTogether(:,3).^2)) ;                    
                                   
%% ========================================================================
%% Fix optimal parameters to get this.  
% Plotting the correlation function
        ZplotCorr = [score2(1:10); score2(11:20); score2(21:30); score2(31:40);  score2(41:50);  score2(51:60)]; 
        score300 = sqrt(sum(((ExptData3)).^2)); 
% Now calculate the R-factor using the best-fit. 
for jj = 1:size(stepnr,2)
        Rfactor(jj) = sqrt(sum((BestScale*dProtMemVacInterp(:,jj).*Contrast - ExptData3).^2))/score300*100;
end
        Zplot2 = [Rfactor(1:10); Rfactor(11:20); Rfactor(21:30); Rfactor(31:40);  Rfactor(41:50); Rfactor(51:60)]; 
% Print to the screen an indication of the quality of the fit. si
        [min(Rfactor) min(scoreComp1) max(score2)]

%% Uncertainties ======================================
close all
% State 2
% plot(sum(Zplot2(:,6:8),2)/3)
% hold on
% plot(sum(Zplot2(3:6,:),1)/3)

% State 1
plot(sum(Zplot2(:,3:4),2)/3)
hold on
plot(sum(Zplot2(2:5,:),1)/3)

%% ===================================================================
% Save theoretical curves for changes to Figures 6C & 6D
        
% BestState1AllFactors = (Best3(BestIndex,:)); 
% BestState1NoMicelle = (Best3(BestIndex,:));         % R-factor 49
 % BestState1NoMicelleNoCrossTerm = (Best3(BestIndex,:));  %       R-factor 52.5  
 % BestState2AllFactors = (Best3(BestIndex,:)); 
 % BestState2NoMicelle = (Best3(BestIndex,:));         % R-factor 25.9
% BestState2NoMicelleNoCrossTerm = (Best3(BestIndex,:));   %      R-factor 29.6
 
 
%% Plot the results of this analysis. 
figure('Position', [750 50 500 400])
contour(Zplot2,[0:5:100],'Linewidth',3)
% contour(Zplot2,[0:10:200],'Linewidth',3)
xlabel('\gamma (Helix E & F)')
ylabel('\delta (Helix C)')
xticks([1:1:10]); xticklabels({'0','1/6','1/3','1/2','2/3','5/6','1','7/6','4/3','1 1/2'})
yticks([1:1:6]); yticklabels({'0','1','2','3','4','5'})
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
colorbar

figure('Position', [350 50 500 400])
plot(q1(FitStart:FitEnd),q1(FitStart:FitEnd).*transpose(ExptData3),'o','LineWidth',0.5,'MarkerSize',4,'color','k');
hold on; 
% plot(q1(FitStart:FitEnd),q1(FitStart:FitEnd).*BestState1NoMicelleNoCrossTerm,'LineWidth',4,'color',[0.2500 0.250 0.95]);
plot(q1(FitStart:FitEnd),q1(FitStart:FitEnd).*(Best3(BestIndex,:)),'LineWidth',4,'color',[0.8500 0.3250 0.0980]);
axis([0.0 1 -6*10^-4 6*10^-4]);
xlabel('q (Å^{-1})')
ylabel('q \cdot \DeltaS(q)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
xticks([0:0.2:1]); 

%% ========================================================================
% Functions used to load difference data from simulations. 
function DataTemp = LoadDiffData(yi)
run_name = yi; 
excited.protmem = load(char(run_name(2)+'_intensities_protmem.mat'));
excited.membrane = load(char(run_name(2)+'_intensities_membrane.mat'));
excited.protein = load(char(run_name(2)+'_intensities_protein.mat'));
ground.protmem = load(char(run_name(1)+'_intensities_protmem.mat'));
ground.membrane = load(char(run_name(1)+'_intensities_membrane.mat'));
ground.protein = load(char(run_name(1)+'_intensities_protein.mat'));
DI(1,:,:) = excited.protmem.intSol(:,1:496) - ground.protmem.intSol(:,1:496);       
DI(2,:,:)  = excited.protmem.intVac(:,1:496) - ground.protmem.intVac(:,1:496);       
DI(3,:,:) = excited.membrane.intSol(:,1:496) - ground.membrane.intSol(:,1:496);   
DI(4,:,:)  = excited.membrane.intVac(:,1:496) - ground.membrane.intVac(:,1:496);    
DI(5,:,:)  = excited.protein.intSol(:,1:496) - ground.protein.intSol(:,1:496);      
DI(6,:,:)  = excited.protein.intVac(:,1:496) - ground.protein.intVac(:,1:496);
DataTemp = DI; 
end
function U17conv = ConvU17(yi,yii,yiii) 
       U17 = yiii; 
       scaleU17 = sum(U17(:,2)); 
       U17CoM = sum(U17(:,1).*U17(:,2))/scaleU17;  
       Econv = 1+(U17(:,1)-U17CoM)/U17CoM; 
       q = yii; 
      for i = 1:round(size(U17,1)/5,0)-1
       TMP(i,:) = U17(5*i,2)*interp1(q*Econv(5*i,1),yi,q(1:135));
      end
      U17conv = [sum(TMP(:,:))/scaleU17  zeros(1,66)];      
%      close all
%      plot(q(1:135),U17conv*max(dProtVac(24,:))/max(U17conv),'color','k')
%      hold on
%      plot(q,dProtVac(24,:),'color','r')
end 