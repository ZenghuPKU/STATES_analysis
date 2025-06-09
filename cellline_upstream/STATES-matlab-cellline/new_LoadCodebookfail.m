function [ geneToSeq, seqToGene ] = new_LoadCodebook( inputPath, remove_index, doReverse )
% new_LoadCodebook

    % load file
    fname = fullfile(inputPath, 'genes.csv');
    f = readmatrix(fname, 'OutputType', 'string', "Delimiter", ',');

    % load gene name and sequence 
    % f(:,1) - gene name, f(:,2) - gene barcode 
    if doReverse
        f(:,2) = reverse(f(:,2));
    end

    for i=1:size(f, 1)
        f(i,2) = new_EncodeBases(f(i,2));
    end

    if ~isempty(remove_index)
        % Sort remove_index in descending order
        remove_index = sort(remove_index, 'descend');

        % Remove characters at specified indices
        for j = 1:length(remove_index)
            f(:,2) = eraseBetween(f(:,2), remove_index(i), remove_index(i));
        end
        segments = {};
        previous_index = 1;
        for i = 1:length(remove_index)
            segments{end+1} = extractBefore(f(:,2), remove_index(i) - (i-1));  % i-1 是因为之前删除操作导致索引位移
            f(:,2) = extractAfter(f(:,2), remove_index(i) - i);
        end
        segments{end+1} = f(:,2);  % 最后一个段
    
        % 反转段的顺序并连接
        f(:,2) = strjoin(flip(segments), '');

    end

    seqToGene = containers.Map(f(:,2), f(:,1));
    geneToSeq = containers.Map(f(:,1), f(:,2));
end
