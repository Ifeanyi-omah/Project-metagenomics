for d in Bat4L Chicken4L control4L Dog4L Egret4L Goat4L Lizard4L Pig4L Rodents4L Sheep4L Squirrel4L
do 
echo $d
nextflow run main.nf --run_dir /localdisk/home/shared/$d --outdir /localdisk/home/shared/nigeria_metagenomics/scylla_output/$d/ -profile docker -resume --local -w work4L > run_$d.log
nextflow run main.nf --run_dir /localdisk/home/shared/$d --outdir /localdisk/home/shared/nigeria_metagenomics/scylla_output/$d/ -profile docker -resume --local -w work4L >> run_$d.log
tar czf /localdisk/home/shared/nigeria_metagenomics/scylla_output/$d.tar.gz $d
rm -rf /localdisk/home/shared/nigeria_metagenomics/scylla_output/$d
rm -rf work4L
done
