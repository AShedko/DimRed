using NIfTI, ExcelReaders, DataFrames
using Revise
using KurchatovFMRI, ProgressMeter, OnlineStats
include("consts.jl")

function main(patient)# Regex gives        year         -  month    -    day
  kd = KData(PATH,"$patient/swvugaf[0-9][0-9][0-9][0-9]\-[0-9][0-9]\-[0-9][0-9].nii","logs/$patient.*", x->x)
  len = size(kd.data,4)
  Ys = reshape(kd.data,(SIZE,len))
  Ss = [Series(LinReg(1)) for i in 1:SIZE] #Sequences
  n_chunks = 600
  low = 1;
  chunksize = len÷n_chunks;
  high = low + chunksize-1;
  Xs = Array{Int}((len,1))
  Xs[:,1] = 1:len

  low = 1;
  high = low + chunksize-1;
  @showprogress 3 "Computing chunks " for ll in 1:n_chunks
    @inbounds for ind in 1:SIZE
      fit!(Ss[ind],( view(low:high, :, :),[Ys[ind,l] for l in low:high]))
    end
    low = low+chunksize
    high = high+chunksize;
    if ll == n_chunks
      high = len
    end
  end

  assert(mean((o->value(o)[1]).(Ss)) > 1e-10, "Mean trend should be significant")
  info("writing trend")
  open("out/trend_$patient.ser", "w") do f
    m = [value(o)[1][1] for o in Ss]
    serialize(f, m)
  end
end

main(parse(Int,ARGS[1]))
#
# matched_data = glob("$patient/swvugaf[0-9]*.nii",PATH)
# data = map(x->niread(joinpath(PATH,x),mmap=true),matched_data)
