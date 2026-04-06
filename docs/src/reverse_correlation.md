# Reverse Correlation

Probabily the most important analysis, the foundation of estimating the spatiotemporal receptive fields, of the first-order component ([Spike-Triggered Average](@ref)) and higher-order component ([Spike-Triggered Covariance](@ref))[^schwartzSpiketriggeredNeuralCharacterization2006].
The basic idea is the features of the stimulus that drive the cellular responses (e.g., spikes) is the features of the receptive fields.

Usually the features here are continous random vairables
(e.g., noise stimuli), it can also be the movie stimulus itself. No matter what feature you choose, the stimulus matrix for this analysis should be organized as `[nDimensions x nTimepoints]`; and the reverse correlation would be `[nDimensions x nFrames]`.

## Spike-Triggered Average

The first order component, which is equivalent to the weighted average of the stimulus before spikes.

```math
\operatorname{STA}(n) = \frac{1}{N} \sum_{i=1}^N \vec{s}(n, t_i)
```

where $t_i$ is the spike time and $\vec{s}(n, t_i)$ is the stimulus feature vector related to the spike.

### bootstrapped spike-tiriggered average

TODO

## Spike-Triggered Covariance

TODO

## References

[^schwartzSpiketriggeredNeuralCharacterization2006]: Schwartz, O., Pillow, J. W., Rust, N. C. & Simoncelli, E. P. Spike-triggered neural characterization. Journal of Vision 6, 484–507 (2006).
