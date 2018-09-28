%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CECI EST LE CODE PRINCIPAL POUR CALCULER TE À PARTIR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  D'UN SPECTRE D'UN PLASMA D'HELIUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
%% ---------- Notes d'utilisation ----------
%====== Les sous fonctions utilisées sont: ======
%   TeTauxDeReactions
%   IntensitePicEdMax.m
%   TeIntTheo_SANS_Input.m
%   TeInitialisation
%   TeEscapeFactor.m 
%   TeDepopulation.m
%   TePopulation.m 
%   TeGraphes1.m
%   TeGraphes2.m

%====== Les variables préprogrammées utilisées sont: ======
%   En2px_2py.mat
%   Allrate1s_2p.mat 
%   AllrateGround_2p.mat
%   AllIntegral_new.mat

%Tous ces fichiers ainsi que les sections efficaces sont dans le sous-dossier "Complements"
%------------------------------------------
%%  ###########################################################################################
%%  ######################################### INPUTS ##########################################
%%  ###########################################################################################

% Range des températures électroniques scannées (eV) (pas plus bas que 0.07 eV)
    Te =0.07:0.01:3;   
    Tg = 300;    
    longueur = 0.01; %longueur auto absorp

%% ###########################################################################################
%% ###################################### DÉBUT DU CODE ######################################
%% ###########################################################################################
 
%% ============== Préchargement de variables ==============
   %Valeurs de beta et de l'intégrale de convolution trouvé à l'aide de TeBeta.m
   load AllIntegral_new.mat
   %Taux de réactions pour Te=0.05eV jusquà 10eV par step de 0.001eV
   load AllratesmetaCroise.mat
   load AllratesmetaResonant.mat
   load Allratesmeta.mat
   load AllratesGround.mat  

addpath('C:\Users\Florence\Dropbox\resultats\no-sample\Calib-He-pur\raies-maxetintegrer')

Intensite=dlmread('ratios1100ppm.txt'); % intensité lumineuse
n23S=dlmread('Meta_1100ppm.txt'); % densité de métastables

Temps=[1:length(Intensite(:,1))];
SingTheo=zeros(length(Temps),length(Te));
TripTheo=zeros(length(Temps),length(Te));

  P=760*133.322;  % Pression

%% ============== Raies d'émission utilisées dans le calcul ==============
    Raies=[388.86 501.57 587.56 667.82 706.52 728.13]'; 
 
% ============== Boucle d'initialisation sur Te pour la sélection des bons taux de réaction ==============
    %Préallocation d'espace
    rateGround_excite=zeros(6,length(Te));          
    rateMet_excite=zeros(4,6,length(Te));
    
for j=1:length(Te)       %Boucle sur Te
    %D'abord on trouve la position du Te correspondant dans Allrates
    position=round((Te(j)/0.001)-49);    %-49 pour retomber à une position 1 pour Te=0.05 et matcher le pas de TeTauxdeReaction
    
    %Ensuite on sélectionne les taux du fondamental aux excités (n>2)
    for l=1:6
        rateGround_excite(l,j)=AllratesGround(l,position);
    end 
    
    %Finalement ceux des n=2 aux n>2
    for k=1:4                       %Boucle sur les n=2
        for l=1:6                   %Boucle sur les n>2
            rateMet_excite(k,l,j)=AllratesmetaCroise(k,l,position);
        end
    end
end

 clear position
 %disp('Taux de réaction selectionnés')

%% ########################### DÉBUT DE L'ANALYSE THÉORIQUE DU CODE ##########################
%% Préallocation d'espace
I_theo      =zeros(length(Temps),length(Te),length(Raies));
Thetaij     =zeros(length(Te),length(Raies));
Doppler     =zeros(length(Te),length(Raies));
VanDerWaals =zeros(length(Te),length(Raies));
Resonant    =zeros(length(Te),length(Raies));

GainFond   =zeros(length(Te),6);
Gain1s     =zeros(length(Te),6);
GainColl2p =zeros(length(Te),6);
GainRadTrap=zeros(length(Te),6);
PerteRad   =zeros(length(Te),6);
PerteColl2p=zeros(length(Te),6);
PerteColl1s=zeros(length(Te),6);
densite    =zeros(length(Te),6);
ContributionFond    =zeros(length(Te),6);

mecanismes=zeros(length(Te)); 
    
for k=1:length(Temps)

   % Calcul du bilan de population: gains=pertes
   [density,Gains,Pertes,Emission,PopFond,Mecanisms]=TeIntTheo_He(n23S(k),0,0,Te,rateGround_excite,rateMet_excite,Tg,longueur,AllIntegral,P,1);

   %% Extraction et mise en mémoire des données pour chaque fichier, densité de metastable et Te
   for i=1:length(Te)
        for m=1:length(Raies)                       % Pour chaque raie...
            I_theo(k,i,m)=Emission(1,i,m);          % Intensité théorique           
        end   
   end
   %Fin Boucle Nm

end
I_theo=abs(I_theo);
%%intensités expérimentales%%

Sing=Intensite(:,1);
Trip=Intensite(:,2);

for i=1:length(Temps)  % on vérifie que les données sont bonnes  
    if Sing(i)>2.2
        Sing(i)=0;
    end
    if Trip(i)>0.15
        Trip(i)=0;
    end
end

TripSingTheo=zeros(length(Temps),length(Te));
TripSing=zeros(length(Temps),length(Te));
TripExpTheo=zeros(length(Temps),length(Te));
SingExpTheo=zeros(length(Temps),length(Te));

Sing1=zeros(length(Temps),length(Te));
Trip1=zeros(length(Temps),length(Te));

for i=1:length(Te)
    for j=1:length(Temps)  
    SingTheo(j,i) = I_theo(j,i,4)./I_theo(j,i,6);
    TripTheo(j,i) = I_theo(j,i,3)./I_theo(j,i,5);
    TripExpTheo(j,i)=Trip(j)./TripTheo(j,i);
    SingExpTheo(j,i)=Sing(j)./SingTheo(j,i);    
    TripSingTheo(j,i) = TripTheo(j,i)/SingTheo(j,i);
    TripSing(j,i)=Trip(j)./Sing(j);
    Sing1(j,i)=Sing(j);
    Trip1(j,i)=Trip(j);
    end
end 

LigneTheo = ones(length(Te),1);

figure   %visualisation des résultats
hold on
yyaxis left
plot(Te,LigneTheo,'-k') 
plot(Te,TripExpTheo(6,:),'-')
ylabel('Ratio Exp/Ratio Théo Triplets')
axis([0 3 0 2])
yyaxis right
hold on
plot(Te,SingExpTheo(6,:),'-')
axis([0 3 0 2])
xlabel('Te') 
ylabel('Ratio Exp/Ratio Théo Singulets')
hold off

TeTripSing=zeros(length(Temps),1);
TeSing=zeros(length(Temps),1);
TeTrip=zeros(length(Temps),1);


for i=1:length(Temps)    % température électronique
   [minimum,position]=min(abs(1-SingExpTheo(i,:)));  %singulets
   if (i>=2) && (ismember(0,SingExpTheo(i,:))==0)              
   k=0;
   while (abs(Te(position)-TeSing(i-1)) > 0.1) & (k<10000)
   SingExpTheo(i,position)=0.1;  
   [minimum,position]=min(abs(1-SingExpTheo(i,:)));  
   k=k+1;
   end
    else 
       TeSing(i)=0;
   end
   TeSing(i)=Te(position);
  
      
   [minimumt,positiont]=min(abs(1-TripExpTheo(i,:)));  %triplets
   if (i>=2) && (ismember(0,TripExpTheo(i,:))==0)
   k=0;
   while (abs(Te(positiont)-TeTrip(i-1)) > 0.1)  & (k<10000)
   TripExpTheo(i,positiont)=0.1;  
   [minimum,positiont]=min(abs(1-TripExpTheo(i,:)));  
    k=k+1;
   end  
   else
       TeTrip(i)=0;
   end
    TeTrip(i)=Te(positiont);   
end

DeltaTe=TeSing-TeTrip;

%%
figure   % Te finale
plot(Temps,TeSing,'ro-')
hold on
plot(Temps,TeTrip,'bd-')
%axis([100 800 0.3 0.6])
legend('Te Sing','Te Trip')
xlabel('Frame')
ylabel('Te')

clear NumFichier i j k l m fileID formatSpec fileID2 formatSpec2 IntExp TeGraph2 TeGraph1 AllIntegral TeGrapheFinaux nbrderaie
