function obj = new_FilterReads_Duo( obj, endBases, showPlots )
% FilterReads_Duo

    % Filter reads by whether they are in the codebook
    % This is just as a sanity check; reads are actually filtered
    % by whether they are in the codebook
    filtBases = obj.basecsMat;
    correctSeqs = zeros(size(filtBases,1), 2); % count sequences that are of correct form

    cs_front = filtBases(:, 1:5);
    cs_back = filtBases(:, 6:9);

    bases_front = new_DecodeCS(cs_front, endBases(1));
    bases_back = new_DecodeCS(cs_back, endBases(2));
    
    disp('bases_front:');
    disp(bases_front);
    disp('bases_back:');
    disp(bases_back);
    
    for i=1:numel(bases_front)
        currSeq = bases_front{i};
	if currSeq(1) == endBases(1) && currSeq(end) == endBases(3)
            correctSeqs(i, 1) = 1;
        end
    end

    for i=1:numel(bases_back)
        currSeq = bases_back{i};
	if currSeq(1) == endBases(2) && currSeq(end) == endBases(4)
            correctSeqs(i, 2) = 1;
        end
    end

    correctSeqs = correctSeqs(:,1) .* correctSeqs(:,2);

    fprintf('Filtration Statistics:\n');
    % fprintf(obj.log, 'Filtration Statistics:\n');
    score_1 = sum(correctSeqs)/size(filtBases,1);
    s1 = sprintf('%f [%d / %d] percent of good reads match barcode pattern %sNNNN%s - %sNNN%s\n',...
        sum(correctSeqs)/size(filtBases,1),...
        sum(correctSeqs),...
        size(filtBases,1),...
        endBases(3),...
        endBases(1),...
	endBases(4),...
	endBases(2));
    fprintf(s1);

    % filter reads based on codebook
    Nreads = numel(obj.allReads);
    inCodebook = zeros(Nreads, 1);
    codebookSeqs = obj.barcodeSeqs;

    for s=1:Nreads
        str = obj.allReads{s};
        %disp(['Read Sequence ', num2str(s), ': ', str]);
        %disp(['Codebook Sequence ', codebookSeqs]);
        if ismember(str, codebookSeqs)         
            inCodebook(s) = 1;
        else
            inCodebook(s) = 0;
        end
    end

    % sum(inCodebook)

    readsToKeep = inCodebook==1;
    obj.goodSpots = obj.allSpots(readsToKeep,:);
    obj.goodReads = obj.allReads(readsToKeep);
    obj.goodScores = obj.allScores(readsToKeep,:);

    if showPlots
        figure(1);
        errorbar(mean(obj.goodScores), std(obj.goodScores),'ko-'); 
        xlim([0 obj.Nround+1]); 
        xlabel('Round'); ylabel('Average qual score');
    end

    score_2 = sum(readsToKeep)/Nreads;
    s2 = sprintf('%f [%d / %d] percent of good reads are in codebook\n',...
        sum(readsToKeep)/Nreads,...
        sum(readsToKeep),...
        Nreads);
    fprintf(s2);

    score_3 = sum(readsToKeep)/sum(correctSeqs);
    s3 = sprintf('%f [%d / %d] percent of form matched reads are in codebook\n',...
        sum(readsToKeep)/sum(correctSeqs),...
        sum(readsToKeep), ...
        sum(correctSeqs));            
    fprintf(s3);

    obj.FilterScores = [score_1 score_2 score_3]; 
    
    if ~isempty(obj.log)
        fprintf(obj.log, s1);
        fprintf(obj.log, s2);
        fprintf(obj.log, s3);
    end
        
end
