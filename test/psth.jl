@test length(spike_histogram(Any[], 0:0.1:1)) == 10
@test sum(spike_histogram(Any[], 0:0.1:1)) == 0
@test sum(spike_histogram(Any[0.5], 0:0.1:1)) == 1
@test sum(spike_histogram(Float32[0.5], 0:0.1:1)) == 1
@test eltype(spike_histogram(Float32[0.5], 0:0.1:1)) <: Float32
@test sum(spike_histogram([Float32[0.5], Float32[0.8]], 0:0.1:1)) == 1
