classdef AdaptiveAmk

    properties (SetAccess = protected)
        trainX
        trainY
        standardizedX
        standardizeStats
        options
        bestK
        bestGCV
    end

    methods (Access = public, Static = false)

        function obj = AdaptiveAmk(trainX, trainY, options, bestK, bestGCV)
            defaultOptions = struct();
            defaultOptions.standardize = true;
            defaultOptions.rangeK = [5,500];
            defaultOptions.epsilon = 1e-6;
            defaultOptions.fastComputation = true;
            defaultOptions.subsetSize = 2000;
            defaultOptions.predBlockSize = 1000;
            defaultOptions.gssOptions = struct(tol=1, extendBoundary=[false, true], ftol=1e-3, verbose=false);

            if exist('options','var')
                obj.options = options;
                optionNames = fieldnames(options);
                if ~ismember("standardize",optionNames)
                    obj.options.standardize = defaultOptions.standardize;
                end
                if ~ismember("rangeK",optionNames)
                    obj.options.rangeK = defaultOptions.rangeK;
                end
                if ~ismember("epsilon",optionNames)
                    obj.options.epsilon = defaultOptions.epsilon;
                end
                if ~ismember("fastComputation", optionNames)
                    obj.options.fastComputation = defaultOptions.fastComputation;
                end
                if ~ismember("subsetSize", optionNames)
                    obj.options.subsetSize = defaultOptions.subsetSize;
                end
                if ~ismember("predBlockSize", optionNames)
                    obj.options.predBlockSize = defaultOptions.predBlockSize;
                end
                if ~ismember("gssOptions", optionNames)
                    obj.options.gssOptions = defaultOptions.gssOptions;
                end
            else
                obj.options = defaultOptions;
            end
            if ~exist('bestK','var')
                obj.bestK = [];
            else
                obj.bestK = bestK;
            end

            if ~exist('bestGCV','var')
                obj.bestGCV = [];
            else
                obj.bestGCV = bestGCV;
            end
            obj.trainX = trainX;
            obj.trainY = trainY;
            if obj.options.standardize
                [obj.standardizedX, obj.standardizeStats] = standardizeData(trainX, "self");
            end
            if isempty(obj.bestK)
                [obj.bestK, obj.bestGCV] = computeBestK(obj);
            end

            return
        end

        function [bestK, bestGCV] = computeBestK(obj)
            if obj.options.fastComputation
                if size(obj.trainX, 1) >= obj.options.subsetSize
                    sampleIdx = datasample(1:size(obj.trainX, 1), obj.options.subsetSize);
                else
                    sampleIdx = 1:size(obj.trainX, 1);
                end
            else
                sampleIdx = 1:size(obj.trainX, 1);
            end
            if obj.options.standardize
                X = obj.standardizedX(sampleIdx, :);
            else
                X = obj.trainX(sampleIdx, :);
            end
            y = obj.trainY(sampleIdx);
            fun = @(k)AdaptiveAmk.computeGCV(X, y, k, obj.options.epsilon, obj.options.predBlockSize);
            gssOptions = obj.options.gssOptions;
            out = goldenSearchInteger(fun,obj.options.rangeK(1),obj.options.rangeK(2),gssOptions);
            bestK = out.sol;
            bestGCV = out.val;
            return
        end



        function predY = predict(obj, testX)
            if obj.options.standardize
                testX = standardizeData(testX, obj.standardizeStats);
                X = obj.standardizedX;
            else
                X = obj.trainX;
            end
            predY = AdaptiveAmk.predictInternal(X, obj.trainY, testX, obj.bestK, obj.options.epsilon, obj.options.predBlockSize);
            return
        end
   
    end

   methods (Access = public, Static = true)

       function [predY, tr] = predictInternal(trainX, trainY, testX, bestK, epsilon, predBlockSize)
            predY = zeros(size(testX, 1), 1);    
            nBlocks = ceil(size(testX, 1)/predBlockSize);
            testIdx = 1:min(predBlockSize, size(testX, 1));
            tr = 0;
            for block = 1:nBlocks    
                [~, D] = knnsearch(trainX, testX(testIdx, :), 'K', bestK);
                bandwidth = D(:,bestK)/3 + epsilon;
                weightMatrix = exp(-((pdist2(testX(testIdx, :),trainX)./bandwidth).^2));
                weightMatrix = weightMatrix./sum(weightMatrix,2);
                predY(testIdx) = weightMatrix*trainY;
                tr = tr + sum(diag(weightMatrix, testIdx(1)-1));
                testIdx = (testIdx(end) + 1):min(testIdx(end) + predBlockSize, size(testX, 1));        
            end
            return
        end

        function gcv = computeGCV(trainX, trainY, k, epsilon, predBlockSize)
            [predY, tr] = AdaptiveAmk.predictInternal(trainX, trainY, trainX, k, epsilon, predBlockSize);
            gcv = sqrt(mean(((trainY - predY)/(1-(tr/length(trainY)))).^2));
            return
        end
    end

end
