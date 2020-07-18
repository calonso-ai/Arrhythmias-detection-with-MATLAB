function alarmResult=challenge(recordName,alarm_type)

% Entrada:
%   recordName
%       Archivo a procesar
%   alarmType
%       Cadena en la que se especifica el tipo de arritmia
%             Bradycardia, Tachycardia,
%             Ventricular_Tachycardia, Atrial_Fibrillation
%
%
% Salida:
%   alarmResult
%       Valor entero; '0' si existe una falsa alarma y '1' si el archivo de entrada
%       contiene una arritmia (True alarm)
%          
%
% Written by Linda Eerikäinen April 7, 2015
% Last modification by Carlos A. Alonso, June, 2017


switch alarm_type
    case 'Atrial_Fibrillation'
        loadData
        dataTest
    otherwise


% Selección de señales para la evaluación de caracteristicas
signalArray = signalSelection(recordName);

% Extracción de caracteristicas para la clasifiación
physiological_features = computeFeatures(alarm_type, signalArray);

%Se usan vectores predictivos hallados mediante Random Forest.
switch alarm_type
    case 'Bradycardia'
        load('bradycardiaRF.mat')
        features = [signalArray.sqi_features(1), signalArray.sqi_features(3), ...
            physiological_features];
        class = predict(bradycardiaRF,features);
        alarmResult = str2double(class);
    case 'Tachycardia'
        features = [signalArray.sqi_features(1), physiological_features]; 
        load('tachycardiaRF.mat')
        class = predict(tachycardiaRF,features);
        alarmResult = str2double(class);
    case 'Ventricular_Tachycardia'
        features = [signalArray.sqi_features(1:5), physiological_features];
        load('ventricularTachycardiaRF.mat')
        class = predict(ventricularTachycardiaRF, features);
        alarmResult = str2double(class);

end

end

end


