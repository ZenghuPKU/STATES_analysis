function colorSeq = new_EncodeBases( seq )
% new_EncodeSOLID

    % construct hash table for encoding
    k = {'AT','GT','TT',...
        'AG', 'GG', 'TG',...
        'AA', 'GA', 'TA'};
    v = {3,1,2,1,2,3,2,3,1};
    coding = containers.Map(k,v);
    
    start = 1;
    back = start + 1;
    colorSeq = "";
    while back <= strlength(seq)
        curr_str = extractBetween(seq, start, back);
        colorSeq = colorSeq + coding(curr_str);
        start = start + 1;
        back = start + 1;
    end
    
end

