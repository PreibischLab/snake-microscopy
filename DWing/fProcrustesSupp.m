% [SupDissim,SupData, SupTform] =  fProcrustesSupp(Model, Data, reflectionFlag)
% Data cell of landmark coordinate table
% Model coordinate table



function [SupDissim,SupData, SupTform] =  fProcrustesSupp(Data, reflectionFlag, Model)

if (~exist('reflectionFlag', 'var')) || isempty(reflectionFlag)
    reflectionFlag = true;
end
if (~exist('Model', 'var')) || isempty(Model)
    Model = Data{1,1};
end
    d = Inf;
    SupDissim = zeros(numel(Data),1);
    SupData = cell(size(Data));
    SupTform = SupData;
    while d > 0.1
        for i=1:numel(Data)
            [SupDissim(i),SupData{i},SupTform{i}] = procrustes(Model,Data{i},'reflection',reflectionFlag);
        end
        
        ModelTmp = sum(cell2mat(reshape(SupData,1,1,[])),3)/numel(Data);
        
        d = mean(sum((ModelTmp-Model).^2,2));
        Model = ModelTmp;
    end

end