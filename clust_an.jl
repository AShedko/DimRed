
push!(LOAD_PATH, "../fMRI_ISC_classify/Jthing")
using ImageView, GtkReactive, Colors
using KurchatovFMRI,ExcelReaders, Vis

mat = open("DimRed/dist_1.ser", "r") do f
  deserialize(f)
end
maxs = [maximum(mat[:,i]) for i in 1:SIZE]
maxs = reshape(maxs,SHAPE);Z = zeros(SHAPE)
kd_m = KData(PATH,"1/mean.nii","logs/1.*",
  x->readxl(DataFrame,x,"Лист1!A1:G$(NEXP+1)"))
mean = kd_m.data
diff = maxs./mean
diff = (x->isnan(x) || x < 0.5  ? 0.0: x).(diff)
diff2 = (x->x>1. ? 0.0 :  x).(diff)
vis(diff2,Z,0)
