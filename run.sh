
for i in *loom; 
do
   bsub -q rerunnable -J loom python loom_convert_to_adata.py $i
done
