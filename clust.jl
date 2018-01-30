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
pointwise(pr,x,y) = pr+abs(abs(x)-abs(y))

"""
Online Kolmogorov distance
"""
kolmogorov(p,x,y) = max(p, abs(abs(x)-abs(y)))

"""
 get_dist_matrix(data)
 ### GDM
 Get Distance matrix for clustering
 data(4d NII) -> sparse distance matrix

```jldoctest
n=10;R = randn(n,n,n,10);
res = dist_matrix(R,1,kolmogorov);

```
"""
function dist_matrix(data,d,metric = pointwise)
  # distm = spzeros(length(data),length(data))
  distm = zeros((2d+1)^3,SIZE)
  # distm = spzeros(SIZE,SIZE)
  @showprogress 3 "Computing Distance matrix" for frame in gen_frame(data,1:LEN÷6)
    # info(CartesianRange(CartesianIndex((0,0,0))+1,CartesianIndex(size(frame))-1))
    @inbounds for i in CartesianRange(CartesianIndex(1,1,1)+d,CartesianIndex(size(frame))-d)
      # i[2]==55 && info(i) # debug
      # for j in δ(i)
      @inbounds for(ind,j) in enumerate(δ(i))
        # distm[s2i(j),s2i(i)] = metric(distm[s2i(j),s2i(i)],frame[i],frame[j])
        distm[ind,s2i(i)] = metric(distm[ind,s2i(i)],frame[i],frame[j])
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


# dbscan(mat, 5.0, 5)
# using PyImport
# @pyimport sklearn.cluster as clst
# dbsc = clst.DBSCAN(eps=20.0,metric="precomputed",n_jobs=-1)
# dbsc[:fit_predict](mat)
