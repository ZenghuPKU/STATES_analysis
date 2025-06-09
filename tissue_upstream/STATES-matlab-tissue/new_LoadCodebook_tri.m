function [geneToSeq, seqToGene] = new_LoadCodebook_tri(inputPath, remove_index, doReverse)
    % new_LoadCodebook

    % Load file
    fname = fullfile(inputPath, 'genes.csv');
    f = readmatrix(fname, 'OutputType', 'string', "Delimiter", ',');

    % Load gene name and sequence
    % f(:,1) - gene name, f(:,2) - gene barcode
    if doReverse
        f(:,2) = reverse(f(:,2));
    end

    for i = 1:size(f, 1)
        f(i, 2) = new_EncodeBases_tri(f(i, 2));
    end
    
    if ~isempty(remove_index)
        % Sort remove_index in ascending order
        remove_index = sort(remove_index, 'ascend');

        % Loop over each row in f
        for row = 1:size(f, 1)
            % Initialize segments for the current row
            segments = strings(1, 0);  % Ensure segments is a string array

            current_sequence = f(row, 2);

            % Split the string into segments based on the adjusted indices
            start_idx = 1;
            for j = 1:length(remove_index)
                end_idx = remove_index(j) - 1;
                if start_idx <= end_idx
                    segments(end+1) = extractBetween(current_sequence, start_idx, end_idx);
                end
                start_idx = remove_index(j) + 1;  % Adjust start_idx to skip the removed character
            end

            % Add the remaining part of the sequence as the last segment
            if start_idx <= strlength(current_sequence)
                segments(end+1) = extractAfter(current_sequence, start_idx - 1);
            end

            % Reorder segments: move the second segment to the front
            if numel(segments) > 1
                reordered_segments = [segments(2), segments(1), segments(3:end)];
            else
                reordered_segments = segments;
            end

            % Concatenate reordered segments
            f(row, 2) = strjoin(reordered_segments, '');
        end
    end
    
    disp(f(:, 2));
    % Create the mappings
    seqToGene = containers.Map(f(:, 2), f(:, 1));
    geneToSeq = containers.Map(f(:, 1), f(:, 2));
end
