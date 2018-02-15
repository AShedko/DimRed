push!(LOAD_PATH, "../fMRI_ISC_classify/Jthing")
using ImageView, GtkReactive, Colors,OffsetArrays
using KurchatovFMRI,ExcelReaders, ImageFiltering
include("consts.jl")

function vis(mri_, mask, thresh=0.1,rng = 0.:0.01:1.)
  mri = Gray.(copy(mri_)/ maximum(mri_))
  mriseg = RGBA.(mri)
  sl = slider(rng)
  mriseg[mask .> thresh] += RGBA(0.,0.,1.,0.5)
  zr, slicedata = roi(mri, (1,2))
  gd = imshow_gui((200, 200), slicedata, (1,2))
  imshow(gd["frame"][1,1], gd["canvas"][1,1], mri, nothing, zr, slicedata)
  imshow(gd["frame"][1,2], gd["canvas"][1,2], mriseg, nothing, zr, slicedata)
  showall(gd["window"])
end

diff = open("DimRed/out/dist_1_S.ser", "r") do f
  deserialize(f)
end
# maxs = [maximum(mat[:,i]) for i in 1:SIZE]
# maxs = reshape(maxs,SHAPE);Z = zeros(SHAPE)
kd_m = KData(PATH,"1/mean.nii","logs/1.*",
  x->x)
mn = kd_m.data

deltas = [Kernel.sobel((true,true,true),x)[1] for x in 1:3]
G = zeros(SIZE)
for i =1:SIZE
  D = [sum(del.*OffsetArray(reshape(diff[:,i], (3,3,3)),-1:1,-1:1,-1:1)) for del in deltas]
  G[i] = sqrt(sum((x->x^2).(D)))
end

Gl = reshape(G,SHAPE)

G0 = (x->isnan(x) ? 0.0: x).(Gl./mn)
G0[mn.<1.] .= 0.
# G2 = (x->x>1. ? 0.0 :  x).(G0)

vis2(mn,G0.>0.01,Gl.>0.1,0.0)
