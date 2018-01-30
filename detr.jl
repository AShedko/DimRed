using NIfTI, ExcelReaders, DataFrames
using Revise,JLD, Glob
using KurchatovFMRI, ProgressMeter, OnlineStats
include("consts.jl")

function main(patient)# Regex gives        year         -  month    -    day   _ time...
  kd = KData(PATH,"$patient/swvugaf[0-9][0-9][0-9][0-9]\-[0-9][0-9]\-[0-9][0-9]_[0-9]*.nii","logs/$patient.*",
    x->x)
  Ys = kd.data
  Ss = [LinReg(1) for i in 1:SIZE] #Sequences
  len = size(kd.data,4)
  n_chunks = 600
  low = 1;
  chunksize = len÷n_chunks;
  high = low + chunksize-1;
  Xs = Array{Int}((len,1))
  Xs[:,1] = 1:len
  chunks = [Array{Float64}((chunksize,1)) for i =1:n_chunks-1]
  push!(chunks, Array{Float64}((len % chunksize,1)))
  for ch in chunks
    ch[:,1] = low:high
    low = low+chunksize
    high = min(high + chunksize, len);
  end
  info("Chunks formed")
  
  low = 1;
  high = low + chunksize-1;
  @showprogress 3 "Computing chunks " for ch in chunks
    @inbounds for ind in 1:SIZE
      Series((ch,[Ys[l][ind] for l in low:high]),Ss[ind])
    end
    low = low+chunksize
    high = min(high + chunksize, len);
  end

  info("writing trend")
  open("trend_$patient.ser", "w") do f
    m = (o->o.β[1]).(Ss)
    serialize(f, m)
  end
end

main(parse(Int,ARGS[1]))
#
# matched_data = glob("$patient/swvugaf[0-9]*.nii",PATH)
# data = map(x->niread(joinpath(PATH,x),mmap=true),matched_data)
