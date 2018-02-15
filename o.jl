using NIfTI
using Revise
using KurchatovFMRI
include("consts.jl")

ni = niread(PATH_Fmt(2) *"swvugaf2016-10-28_10-59-114454-03620-03620-1.nii" )
show(length(ni))
dat = ones(Int16,(91, 109, 91))
hdr = ni.header
hdr.datatype = NIfTI.nidatatype(Int16)
res = NIfTI.NIVolume{Int16, 3, Array{Int16,3}}(hdr, ni.extensions, dat)
niwrite("out/test.nii",res)
