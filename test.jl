using Base.Mmap, NIfTI, ResumableFunctions
push!(LOAD_PATH, "./fMRI_ISC_classify/Jthing")
import ISCLib

@resumable function getind(data)
  l = size(data)[4]
  for i=1:l
    @yield data[:,:,:,i]
  end
  [0,0]
end


function get_n(subj::Int, ind::Int)
    fnames = readdir(PATH_Fmt(subj))[1:end] # first is detrended,second elem is mean
    return niread(joinpath(PATH_Fmt(subj),fnames[ind]), mmap = true)
end

ni = get_n(1,1)

function g(ni)
  s=0.0 ;
  for im in getind(ni)
    !isnull(im) ? s+=im[length(im) ÷ 2] : break;
  end
  s
end


using TSne, MNIST

rescale(A, dim::Integer=1) = (A .- mean(A, dim)) ./ max.(std(A, dim), eps())

data, labels = traindata()
data = convert(Matrix{Float64}, data[:, 1:2500])'
# Normalize the data, this should be done if there are large scale differences in the dataset
X = rescale(data, 1)

Y = tsne(X, 2, 50, 1000, 20.0)

using Gadfly
theplot = plot(x=Y[:,1], y=Y[:,2], color=string.(labels[1:size(Y,1)]))
draw(PDF("myplot.pdf", 4inch, 3inch), theplot)
