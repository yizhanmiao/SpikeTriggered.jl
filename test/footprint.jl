
@test SpikeTriggered.collapse_footprint(
    [
        0 0 0; 0 1 0; 0 0 0;;;
        0 2 0; 0 2 0; 0 0 1
    ];
    gridsize=3,
) == [
    0 2 0;
    0 1 0;
    0 0 1
]
