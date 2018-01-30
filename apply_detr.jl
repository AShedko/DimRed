using NIfTI
using Revise,JLD, Glob
using KurchatovFMRI, ProgressMeter, OnlineStats
include("consts.jl")

function main(patient)
  kd_original = KData(PATH,"1/detre*.nii","logs/1.*",
    x->x)
  ni = kd_original.data
  kd = KData(PATH,"$patient/swvugaf[0-9][0-9][0-9][0-9]\-[0-9][0-9]\-[0-9][0-9]_[0-9]*.nii","logs/$patient.*",
    x->x))
  Ys = kd.data
  len = size(kd.data,4)
  βs = open("trend_$patient.ser", "r") do f
      deserialize(f, m)
  end
  res = zeros(SIZE, len)
  for i = 1 : len
    for v in eachindex(Ys[i])
      res[v,i] = Ys[i][v]-βs[v]
    end
  end

  res = NIfTI.NIVolume(ni.header, ni.extensions, reshape(res, (SHAPE..., len)) # Cheap, as reshape returns a view
  niwrite("out/$pref|segmentation_result.nii",res)
end

main(parse(Int,ARGS[1]))
