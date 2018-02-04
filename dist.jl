using NIfTI, ExcelReaders, DataFrames
using Revise,JLD
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
      @inbounds for(ind,j) in enumerate(δ(i))
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
function dist_matrix_onl(data,d,metric = pointwise)
  # distm = spzeros(length(data),length(data))
  distm = zeros((2d+1)^3,SIZE)
  # exponentially weighted Mean
  series = [Series(ExponentialWeight(0.3), Mean())  for i in 1:SIZE]
  for (i,el) in enumerate(data[:,:,:,1])
    fit!(series[i],el) # initial values for series
  end
  @showprogress 3 "Computing Distance matrix" for frame in gen_frame(data,2:LEN)
    @inbounds for i in CartesianRange(CartesianIndex(1,1,1)+d, CartesianIndex(size(frame))-d)
      fit!(series[i],frame[i]); # fill online stat for current point
      ii = s2i(i);
      @inbounds for(ind,j) in enumerate(δ(i))
        distm[ind,ii] = metric(distm[ind,ii],series[ii],series[ii])
      end
    end
  end
  distm
end

function main(patient)
  kd = KData(PATH,"$patient/detre*.nii","logs/$patient.*",
    x->readxl(DataFrame,x,"Лист1!A1:G$(NEXP+1)"))
  open("dist_$patient.ser", "w") do f
    @time m = dist_matrix(kd.data,2,kolmogorov)
    serialize(f, m)
  end
end

main(parse(Int,ARGS[1]))
