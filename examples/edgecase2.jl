using LaTeXStrings
using OptimalTrainControl
using Plots

trackX = [0,2e3,3e3,5e3]
trackY = [0,0,60,65]

track = HillyTrack(trackX, trackY)

params = ModelParams(track = track)

segs = OptimalTrainControl.getmildsegments(params)

OptimalTrainControl.try_link(segs[3].start,segs[4], :Coast, params)

l = OptimalTrainControl.link(segs[3],segs[4],params)
l[2]