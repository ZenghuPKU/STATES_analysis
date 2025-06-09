function colorSeq = new_EncodeBases_tri( seq )
% new_EncodeSOLID

    % construct hash table for encoding
    k = {'AT','GT','TT',...
        'AG', 'GG', 'TG',...
        'AA', 'GA', 'TA'};
    v = {3,1,2,1,2,3,2,3,1};
    coding = containers.Map(k,v);

    final_k = {'GC','GT'};
    final_v = {1,2};
    final_coding = containers.Map(final_k, final_v);

    start = 1;
    back = start + 1;
    colorSeq = "";

    while back <= strlength(seq)
        curr_str = extractBetween(seq, start, back);
	if back == strlength(seq) % omic
            colorSeq = colorSeq + final_coding(curr_str); 
        else
            colorSeq = colorSeq + coding(curr_str); 
        end
        start = start + 1;
        back = start + 1;
    end
    
end

