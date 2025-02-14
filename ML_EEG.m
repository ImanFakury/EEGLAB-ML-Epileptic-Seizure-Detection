%part1
dataPaths = {
    'C:\Users\imanf\Downloads\WOS1\chb01_02.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_03.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_04.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_15.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_18.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_16.edf',
    'C:\Users\imanf\Downloads\WS1\chb01_26.edf'
};
%converting all data into a matrix
EEGData = cell(1, length(dataPaths));
for i = 1:length(dataPaths)
    EEGData{i} = edfread(dataPaths{i});
end
numericData = cell(1, length(EEGData));
for i = 1:length(EEGData)
    if istimetable(EEGData{i})
        numericData{i} = table2array(EEGData{i});
    else
        numericData{i} = EEGData{i};
    end
end
% converting each second which has 256 elements as a vector because of
%sampeling frequency into 256 main matrix elements
for i = 1:length(numericData)
    data = numericData{i};
    numSamples = size(data, 1) * 256; 
    numChannels = size(data, 2);
    concatenatedData = zeros(numSamples, numChannels);  
    for j = 1:size(data, 1)
        startIdx = (j-1)*256 + 1;
        endIdx = j*256;
        concatenatedData(startIdx:endIdx, :) = cell2mat(data(j, :)); 
    end
    numericData{i} = concatenatedData;
end




%Extracting 14 epochs with 10 minute length accroding to instructions 
%and datasets seizure times
data = numericData{1}; 
intervalLength = 10 * 60 * 256; 
A1 = data(1:intervalLength, :);
A2 = data(intervalLength+1:2*intervalLength, :);
A3 = data(2*intervalLength+1:3*intervalLength, :);
A4 = data(3*intervalLength+1:4*intervalLength, :);
A5 = data(4*intervalLength+1:5*intervalLength, :);
A6 = data(5*intervalLength+1:6*intervalLength, :);
data = numericData{2}; % 'chb01_03.edf'
seizureStart = 2996 * 256; % Seizure start time in samples
seizureEnd = 3036 * 256; % Seizure end time in samples
B1 = data(seizureStart-10*60*256:seizureStart-1, :); 

data = numericData{3}; % 'chb01_04.edf'
seizureStart = 1467 * 256; % Seizure start time in samples
seizureEnd = 1494 * 256; % Seizure end time in samples
B2 = data(seizureEnd:seizureEnd+10*60*256-1, :); 

data = numericData{2}; % 'chb01_03.edf'
C1 = data(seizureStart-5*60*256:seizureStart+5*60*256-1, :); 

data = numericData{3}; % 'chb01_04.edf'
C2 = data(seizureStart-5*60*256:seizureStart+5*60*256-1, :); 

data = numericData{4}; % 'chb01_15.edf'
seizureStart = 1732 * 256; % Seizure start time in samples
seizureEnd = 1772 * 256; % Seizure end time in samples
C3 = data(seizureStart-5*60*256:seizureStart+5*60*256-1, :); 

data = numericData{6}; % 'chb01_16.edf'
seizureStart = 1015 * 256; % Seizure start time in samples
seizureEnd = 1066 * 256; % Seizure end time in samples
C4 = data(seizureStart-5*60*256:seizureStart+5*60*256-1, :); 

data = numericData{5}; % 'chb01_18.edf'
seizureStart = 1720 * 256; % Seizure start time in samples
seizureEnd = 1810 * 256; % Seizure end time in samples
C5 = data(seizureStart-5*60*256:seizureStart+5*60*256-1, :); 

data = numericData{7}; % 'chb01_26.edf'
seizureStart = 1720 * 256; 
seizureEnd = 1810 * 256; 
C6 = data(seizureStart-5*60*256:seizureStart+5*60*256-1, :); 


%part4
%Extracting 16 second epochs(37 epochs)from each one of 14 10 minute
%intervals
matrices = {A1, A2, A3, A4, A5, B1, C1, C2, C3, C4,A6 ,B2,C5, C6};
epochLength = 16 * 256; 
stepSize = epochLength; 
features = cell(1, length(matrices));
for i = 1:length(matrices)
    signal = matrices{i};
    signalLength = size(signal, 1);
    numEpochs = floor(signalLength / stepSize);
    if numEpochs > 37
        numEpochs = 37;
    end
    features{i} = cell(numEpochs, size(signal, 2)); 

    %Extracting 5 features for each channel of 16 second epochs
    for j = 1:numEpochs
        startIdx = (j-1)*stepSize + 1;
        endIdx = startIdx + stepSize - 1;
        if endIdx > signalLength
            endIdx = signalLength;
        end
        epoch = signal(startIdx:endIdx, :);

        if iscell(epoch)
            epoch = cell2mat(epoch);
        end
        for k = 1:size(epoch, 2)
            channelEpoch = epoch(:, k);
            [pxx, freq] = pwelch(channelEpoch, [], [], [], 256);
            psd = 10*log10(pxx);
            psd = psd / sum(psd);
            entropy = -sum(psd .* log2(psd));
            meanVal = mean(channelEpoch);
            stdVal = std(channelEpoch);
            minVal = min(channelEpoch);
            maxVal = max(channelEpoch);
            features{i}{j, k} = [meanVal, stdVal, minVal, maxVal, entropy];
        end
    end
end

%Handeling possible nan values in features matrix
for i = 1:length(features)
    for j = 1:size(features{i}, 1)
        for k = 1:size(features{i}, 2)
            if iscell(features{i}{j, k})
                if all(cellfun(@isnumeric, features{i}{j, k}))
                    features{i}{j, k} = cell2mat(features{i}{j, k});
                else
                    features{i}{j, k} = NaN;
                end
            end
        end
    end
end

%converting features from a 14*23*37*5 matrix into a 
% 14*4255 matrix
vectorizedFeatures = zeros(length(features), 5 * 37 * 23);
for i = 1:length(features)
    featureVector = [];
    for j = 1:size(features{i}, 1)
        for k = 1:size(features{i}, 2)
            if isnumeric(features{i}{j, k})
                featureVector = [featureVector, features{i}{j, k}];
            end
        end
    end
    vectorizedFeatures(i, :) = featureVector;
end



%part 5
%We only use train dataset for feature selection
trainMatrix = vectorizedFeatures(1:6, :);%Non seizure data
testMatrix = vectorizedFeatures(7:10, :);%Seizure data
%Feature extraction using p-value
selectedFeatures = [];
pValues = [];
for i = 1:size(trainMatrix, 2)
    [h, p, ci, stats] = ttest2(trainMatrix(:, i), testMatrix(:, i), 'Alpha', 0.001);
    if p < 0.001
        selectedFeatures = [selectedFeatures, i];
        pValues = [pValues, p];
    end
end

disp('Selected features:');
disp(selectedFeatures);
disp('P-values of selected features:');
disp(pValues);
selectedIndices = selectedFeatures;

%Updating train and test dataset so they have values
%only related to selected features
mySelectedFeatures = zeros(length(matrices), length(selectedIndices));
for i = 1:length(matrices)
    featureVector = vectorizedFeatures(i, :);
    mySelectedFeatures(i, :) = featureVector(selectedIndices);
end


%part 6
%Labeling dataset(-1 for no seizure and 1 for seizure data)
labels = [-1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1];
%We have 10 train data according to project document
X_train = real(mySelectedFeatures(1:10, :)); 
X_test = real(mySelectedFeatures(11:14, :)); 
Y_train = labels(1:10);
Y_test = labels(11:14); 

%I used tic and toc for latency measurment
tic;
%training SVM model by paramethers explaind in document
svmModel = fitcsvm(X_train, Y_train, 'KernelFunction', 'linear', 'Standardize', true);
Y_pred = predict(svmModel, X_test);
svm_latency = toc;
%Computing accuracy of model
numCorrect=0
for i=1:length(Y_test)
    if Y_pred(i)==Y_test(i)
        numCorrect=numCorrect+1;
    end
end
accuracy=numCorrect/4*100;
disp(['Accuracy of SVM classifier: ', num2str(accuracy), '%']);
disp(['Latency of SVM classifier: ', num2str(svm_latency), ' seconds']);


%part 7
%k=2 was best for my data and KNN model 
k = 2 ; 
tic;
%training KNN model by paramethers explaind in document
knnModel = fitcknn(X_train, Y_train, 'NumNeighbors', k, 'Distance', 'minkowski');
Y_pred_knn= predict(knnModel, X_test); 
knn_latency = toc;
%Computing accuracy of model
numCorrect_knn =0;
for i=1:length(Y_test)
    if Y_pred_knn(i)==Y_test(i)
        numCorrect_knn=numCorrect_knn+1;
    end
end
accuracy_knn=numCorrect_knn/4*100;
disp(['Accuracy of KNN classifier: ', num2str(accuracy_knn), '%']);
disp(['Latency of KNN classifier: ', num2str(knn_latency), ' seconds']);

%part 8
truePositives_svm = 0;
falsePositives_svm = 0;
trueNegatives_svm = 0;
falseNegatives_svm = 0;
truePositives_knn = 0;
falsePositives_knn = 0;
trueNegatives_knn = 0;
falseNegatives_knn = 0;

%Iterating in Y-pred and Y-test and Y-pred-knn
%in order to compute number of TP,FM,TN,FP predictins
%in each model
for i = 1:length(Y_test)
    if Y_pred(i) == 1 && Y_test(i) == 1
        truePositives_svm = truePositives_svm + 1;
    elseif Y_pred(i) == 1 && Y_test(i) == -1
        falsePositives_svm = falsePositives_svm + 1;
    elseif Y_pred(i) == -1 && Y_test(i) == -1
        trueNegatives_svm = trueNegatives_svm + 1;
    elseif Y_pred(i) == -1 && Y_test(i) == 1
        falseNegatives_svm = falseNegatives_svm + 1;
    end
    if Y_pred_knn(i) == 1 && Y_test(i) == 1
        truePositives_knn = truePositives_knn + 1;
    elseif Y_pred_knn(i) == 1 && Y_test(i) == -1
        falsePositives_knn = falsePositives_knn + 1;
    elseif Y_pred_knn(i) == -1 && Y_test(i) == -1
        trueNegatives_knn = trueNegatives_knn + 1;
    elseif Y_pred_knn(i) == -1 && Y_test(i) == 1
        falseNegatives_knn = falseNegatives_knn + 1;
    end
end

%Computing specifity and sensivity according to formula
sensitivity_svm = truePositives_svm / (truePositives_svm + falseNegatives_svm);
specificity_svm = trueNegatives_svm / (trueNegatives_svm + falsePositives_svm);
sensitivity_knn = truePositives_knn / (truePositives_knn + falseNegatives_knn);
specificity_knn = trueNegatives_knn / (trueNegatives_knn + falsePositives_knn);
disp(['Sensitivity of SVM classifier: ', num2str(sensitivity_svm)]);
disp(['Specificity of SVM classifier: ', num2str(specificity_svm)]);
disp(['Sensitivity of KNN classifier: ', num2str(sensitivity_knn)]);
disp(['Specificity of KNN classifier: ', num2str(specificity_knn)]);

%Calculating leave one out cross validation
X = real(mySelectedFeatures);
Y = labels;
accuracies_svm = zeros(length(labels), 1);
accuracies_knn = zeros(length(labels), 1);

for i = 1:length(labels)
    X_train = X([1:i-1, i+1:end], :);
    Y_train = labels([1:i-1, i+1:end]);
    X_test = X(i, :);
    Y_test = labels(i);
    svmModel = fitcsvm(X_train, Y_train, 'KernelFunction', 'linear', 'Standardize', true);
    Y_pred_svm = predict(svmModel, X_test);
    accuracies_svm(i) = (Y_pred_svm == Y_test);

    k = 2;
    knnModel = fitcknn(X_train, Y_train, 'NumNeighbors', k, 'Distance', 'minkowski');
    Y_pred_knn = predict(knnModel, X_test);
    accuracies_knn(i) = (Y_pred_knn == Y_test);
end
average_accuracy_svm = mean(accuracies_svm) * 100;
average_accuracy_knn = mean(accuracies_knn) * 100;

disp(['Average accuracy of SVM classifier using LOOCV: ', num2str(average_accuracy_svm), '%']);
disp(['Average accuracy of KNN classifier using LOOCV: ', num2str(average_accuracy_knn), '%']);

