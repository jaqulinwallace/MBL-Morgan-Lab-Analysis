% -------------------------------------------------------------------------
% PCC and MCC Code
% [JW 2023]

% -------------------------------------------------------------------------
% Description:
% Costes Threshold
% Pearson's Correlation Coefficient (PCC) 
% Manders' Coefficient (MCC)
% -------------------------------------------------------------------------

clear all
clc

% User Variables:
species1 = 'pS129 alpha-synuclein'; % Green channel
species2 = 'SV2'; % Red channel
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Code Starts Here:
disp('Batching...')
disp('...')

foldern=uigetdir();
filesn=dir([foldern,'\*.tif']);
remove=zeros(1,numel(filesn));
for i=1:numel(filesn)
    if strfind(filesn(i).name,'Manders')
        remove(i)=1;
    end
end
filesn(logical(remove))=[];

% Load in images:
for j=1:numel(filesn)
    I = tiffreadVolume([foldern,'\',filesn(j).name]);
    disp([num2str(j),'\',num2str(numel(filesn))])
    
    % Split image into channels:
    Ired   = double(I(:,:,1));
    Igreen = double(I(:,:,2));
    
    Red   = reshape(Ired,1,[])';
    Green = reshape(Igreen,1,[])';
    
    % plot(Red,Green,'.')
    
    % Costes Threshold:
    [p,S]=polyfit(Red,Green,1);
    CostesThresh=max(Red); % Start at max threshold value
    Pearson=1;
    tempthresh = zeros(CostesThresh,1)+1;
    temppearson = zeros(CostesThresh,1)+1;

    while Pearson>0
        logicalThresh=Green<(p(1)*CostesThresh+p(2)) & Red<CostesThresh;
        [R,~]=corrcoef(Green(logicalThresh),Red(logicalThresh));
        Pearson=R(1,2);
        temppearson(CostesThresh)=Pearson;            
        if Pearson>0
            tempthresh(CostesThresh)=CostesThresh;
            CostesThresh=CostesThresh-1;
        end                    
    end
    lidxCostes=Green>(p(1)*CostesThresh+p(2)) & Red>CostesThresh;
    
    % Pearson's
    [R,P]=corrcoef(Green(lidxCostes),Red(lidxCostes));
    R_all(j)=R(1,2);
    R(1,2)
    % Manders'
    MandersG(j) = sum(Green(Green>(p(1)*CostesThresh+p(2))))./sum(Green(:));
    MandersR(j) = sum(Red(Red>CostesThresh))./sum(Red(:));
    % -------------------------------------------------------------------------

end

% Final values for batched folder:
R_all(isnan(R_all))=[];
MandersG(isnan(MandersG))=[];
MandersR(isnan(MandersR))=[];

% Figures of PCC and MCC coefficients:
    h1=figure;
    binsn=[-1:0.1:1];
    hist(R_all,binsn)
    title([species1,' vs ',species2,', N=',num2str(numel(R_all)),', mean=',num2str(mean(R_all)),'±',num2str(std(R_all)./sqrt(numel(R_all)))])
    ylabel('Frequency')
    xlabel('Pearson Coefficient')
    set(gca,'Xlim',[-1,1])
    
    h2=figure;
    plot(MandersR,MandersG,'.','MarkerSize',14)
    xlim([0 1])
    ylim([0 1])
    title(['Manders Correlation Coefficients: ',species1,' vs ',species2],'FontSize',14)
    ylabel(['Manders Coefficient: ',species1],'FontSize',14)
    xlabel(['Manders Coefficient: ',species2],'FontSize',14)
    set(gca,'FontSize',14);

% Saving data:
timestamp=datestr(now,'mm-dd-yy+HH-MM-SS');
data = [R_all',MandersG',MandersR'];
saveas(h1,[foldern,'\Pearson_',timestamp])
saveas(h2,[foldern,'\Manders_',timestamp])
save([foldern,'\Manders_PearsonCoefficients_',timestamp],'species1','species2','MandersG','MandersR','CostesThresh','p','R_all','foldern','filesn')
writematrix(data,[foldern,species1, 'MCC-PCCvalues_',timestamp,'.csv'])

disp('Saved.')
disp('...')
disp('Done.')


