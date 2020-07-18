function [tmpTemplate, peakOrientation] = getTemplate(ECG,iR,fs)
    %
    % Forms a QRS template from similar peaks
    %
    % [tmpTemplate, peakOrientation] = getTemplate(ECG,iR,fs)
    %
    % Inputs:
    % ECG = signal
    % iR  = indices of R-peaks
    % fs  = sampling frequency
    %
    % Outputs:
    % tmpTemplate = QRS template
    % peakOrientation  = orientation of the QRS complex in the template
    % 
    % Written by Michiel Rooijakkers, modified by Linda Eerikäinen
    % Last modification by Carlos A. Alonso, May, 2017
    
    QRS    = round(fs * [-0.025 0.025]);
    Delta1 = round(fs * 0.010);
    Delta2 = round(fs * 0.025);
    Delta3 = ceil(fs * 0.050);
    Ni     = 20;
    
    peakIndex  = Delta2 - QRS(1) + 1;
    D1range = peakIndex + (-Delta1:Delta1);
    D2range = peakIndex + (-Delta2:Delta2);
    
    % Initialize template shape based on first Ni Rpeaks using the first Ni
    % peaks which are most similar in shape. 
    
    j = 1;
    
    QRSarray = nan(Ni,diff(QRS)+1+2*Delta2);
    for i = j:Ni
        QRSarray(i,:) = ECG(iR(i)+QRS(1)-Delta2:iR(i)+QRS(2)+Delta2);
        QRSarray(i,:) = QRSarray(i,:) - mean(QRSarray(i,:));
    end
    
    CC = nan(Ni);
    for i = 1:Ni
        CC(i,:) = sum(repmat(QRSarray(i,:),Ni,1) .* QRSarray,2)';
    end
    
    % Remove QRS complexes which deviate from the expected norm a lot.
    AC  = CC(eye(Ni) == 1);
    if std(AC)/quantile(AC,0.95) > 0.12
        selection = find(median(AC) + std(AC) > AC);
    else
        selection = 1:length(AC);
    end
        
    CC  = CC(selection,selection);
    
    % Make a temporary QRS template from all QRS complexes with the same
    % sign (the largest group of complexes with the same sign is selected)
    tmpTemplate = mean(QRSarray(selection(sign(CC(1,:)) == sign(sum(sign(CC(1,:))))),:),1);
   
    
    % Find the exact peak location of the temporary QRS template
    peakOrientation = sign(sum(tmpTemplate(D1range)));
        
    % QRStemplate is initialized as an array of QRS segments which have the
    % same sign, all of which are cut to the correct length around the
    % common newly calculated peak position
    QRSbuffer = nan(length(selection),Delta3);

    for i = 1:length(selection)
        [~,I] = max(peakOrientation*QRSarray(selection(i),D2range));
        peakOffset = I - Delta2 - 1;
        QRSbuffer(i,:) = QRSarray(selection(i),peakIndex + peakOffset + (QRS(1):QRS(2)));
    end
end


