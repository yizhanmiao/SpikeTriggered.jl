@test begin
    spk = [0.1, 0.4, 0.402, 0.5]
    bt = spike_split_burst(spk; detector=spike_detect_burst_lgn)
    bt.tonic == [0.1, 0.5] && bt.burst == [[0.4, 0.402]]
end

@test begin
    spk = [0.1, 0.4, 0.402, 0.5]
    bt = spike_split_cardinal(spk; detector=spike_detect_burst_lgn)
    bt.tonic == [0.1, 0.5] && bt.burst == [0.4]
end

@test begin
    spk = [0.1, 0.4, 0.5]
    bt = spike_split_cardinal(spk; detector=spike_detect_burst_lgn)
    bt.tonic == [0.1, 0.4, 0.5] && bt.burst == []
end