for d in Ms_data failed_run Rodent1a Rodent1b Rodent2a Rodent2b Rodent3a Rodent3b Rodent4 Rodent5
do 
echo $d
nextflow run main.nf --run_dir /localdisk/home/shared/$d --outdir /localdisk/home/shared/nigeria_metagenomics/scylla_output/$d/ -profile docker -resume --local > run_$d.log
nextflow run main.nf --run_dir /localdisk/home/shared/$d --outdir /localdisk/home/shared/nigeria_metagenomics/scylla_output/$d/ -profile docker -resume --local >> run_$d.log
tar czf /localdisk/home/shared/nigeria_metagenomics/scylla_output/$d.tar.gz $d
rm -rf /localdisk/home/shared/nigeria_metagenomics/scylla_output/$d
rm -rf work
done

