using NIfTI, ExcelReaders, DataFrames
using Revise,JLD, OnlineStats
using KurchatovFMRI, ProgressMeter
include("consts.jl")

δ(ind::CartesianIndex, d=1) = CartesianRange(ind-d, ind+d)
i2s(i) = ind2sub(SHAPE,i)
s2i(s) = sub2ind(SHAPE,s...)
s2i(s::CartesianIndex) = sub2ind(SHAPE,s.I...)

@assert i2s(200)|>s2i == 200

"""
Pointwise distance.
Ignores previous observation
"""
pointwise(pr,x,y) = pr + abs(x-y)

"""
Online "Kolmogorov" distance
"""
kolmogorov(p,x,y) = max(p, abs(x-y))

"""
 get_dist_matrix(data)
### GDM
 Get Distance matrix for clustering
 data(4d NII) -> Matrix [ d(vᵢ,δ(vᵢ)) for vᵢ ]

```jldoctest
n=10;R = randn(n,n,n,10);
res = dist_matrix(R,1,kolmogorov);
```
"""
function dist_matrix(data,d,metric = pointwise)
  # distm = spzeros(length(data),length(data))
  distm = zeros((2d+1)^3,SIZE)
  @showprogress 3 "Computing Distance matrix" for frame in gen_frame(data,1:LEN÷6)
    @inbounds for i in CartesianRange(CartesianIndex(1,1,1)+d,CartesianIndex(size(frame))-d)
      @inbounds for (ind,j) in enumerate(δ(i,d))
        distm[ind,s2i(i)] = metric(distm[ind,s2i(i)],frame[i],frame[j])
      end
    end
  end
  distm
end

"""
get_dist_matrix_щтд(data)
### GDM
Get Distance matrix for clustering using OnlineStats with
  EWMA
data(4d NII) -> Matrix [ d(vᵢ,δ(vᵢ)) for vᵢ ]
"""
function dist_matrix_onl(data,d,metric = pointwise, rng = 1:size(data,4))
  # distm = spzeros(length(data),length(data))
  distm = zeros((2d+1)^3,SIZE)
  series = [Series(ExponentialWeight(0.3), Mean()) for i in 1:SIZE]
  # exponentially weighted Mean
  for (i,el) in enumerate(data[:,:,:,1])
    fit!(series[i],el) # initial values for series
  end

  @showprogress 2 "Computing frame" for ll = rng[2:end]
    frame = view(data, :, :, :, ll)
    @inbounds for i in CartesianRange(CartesianIndex(1,1,1)+d, CartesianIndex(size(frame))-d)
      fit!(series[s2i(i)],frame[i]); # fill online stat for current point
      ii = s2i(i);
      @inbounds for (ind,j) in enumerate(δ(i,d))
        distm[ind,ii] = metric(distm[ind,ii],value(series[ii])[1],value(series[s2i(j)])[1])
      end
    end
  end
  return distm
end

function get_range(meta, t::Array{String,1} = ["S1","S2"], rb::Int64=1)
  meta = meta[reduce((x,y)-> x .| y,(Array(meta[:Stimul_1].== w) for w in t)),:] # all rows with stimuls in t
  meta = meta[meta[Symbol("Rest dur")].>0.,:] # No negative lenghts
  ind_end = meta[:Response]
  #ceil returns float64
  dat = vcat(map(x-> collect(Int(round((x-6)*2)) : Int(round((x+2)*2))), ind_end[1:end-1-rb,:] )...)
end

function main(patient, t)
  kd = KData(PATH,"$patient/detre*.nii","logs/$patient.*",
    x->readxl(DataFrame,x,"Лист1!A1:G$(NEXP+1)"))
  open("out/dist_$(patient)_$(t[1][1]).ser", "w") do f
    @time m = dist_matrix_onl(kd.data,1,pointwise, get_range(kd.meta,t))
    serialize(f, m)
  end
end

Stimuli = ["S1","S2","V1","V2","V3","V4"]
if length(ARGS)>=2
  main(parse(Int,ARGS[1]), filter(x->ismatch(Regex(ARGS[2]), x), Stimuli))
elseif  length(ARGS)==1
  main(parse(Int,ARGS[1]), ["S1","S2"])
else
  show("No Args specified! 1st arg - subject_id, 2nd arg - Stimul regex")
end
