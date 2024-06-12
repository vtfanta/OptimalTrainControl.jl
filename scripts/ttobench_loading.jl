using OptimalTrainControl

# Download the track from the TTOBench repository (https://github.com/dkouzoup/TTOBench/tree/main/tracks)
track = load_ttobench_track("CH_Fribourg_Bern.json")

plot(track)