@test begin
    spk = [0.1, 0.4, 0.402, 0.5]
    bt = split_tonic_burst(spk; detector=burst_detect_lgn)
    bt.tonic == [0.1, 0.5] && bt.burst == [[0.4, 0.402]]
end

@test begin
    spk = [0.1, 0.4, 0.402, 0.5]
    bt = split_tonic_cardinal(spk; detector=burst_detect_lgn)
    bt.tonic == [0.1, 0.5] && bt.burst == [0.4]
end

@test begin
    spk = [0.1, 0.4, 0.5]
    bt = split_tonic_cardinal(spk; detector=burst_detect_lgn)
    bt.tonic == [0.1, 0.4, 0.5] && bt.burst == []
end