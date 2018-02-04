using NIfTI
using Revise,JLD, Glob
using KurchatovFMRI, ProgressMeter, OnlineStats
include("consts.jl")

function main(patient)
  kd_original = KData(PATH,"1/detre*.nii","logs/1.*", x->x)
  ni = kd_original.data
  kd = KData(PATH,"$patient/swvugaf[0-9][0-9][0-9][0-9]\-[0-9][0-9]\-[0-9][0-9]\.nii","logs/$patient.*",  x->x)
  len = size(kd.data,4)
  Ys = reshape(kd.data,(SIZE, len))
  βs = open("out/trend_$patient.ser", "r") do f
    deserialize(f)
  end
  res = zeros(Int16,(SIZE, len))
  @showprogress 1 for i = 1 : len
    @inbounds for v =1:SIZE
      res[v,i] = Int16(floor(Ys[v,i]-i*βs[v]))
    end
  end

  # we pass the elemet type to present to the user, the dimension and the backing datastructure
  resn = NIfTI.NIVolume{Int16, 4, Array{Int16,4}}(ni.header, ni.extensions, reshape(res, (SHAPE..., len)))
  niwrite("$PATH/$patient/detrended.nii",resn)
end

main(parse(Int,ARGS[1]))
