config_file=$1

snakemake --snakefile /n/data1/hms/dbmi/park/dominika/testing/smaht/pipelines/MosaicForecast/Snakefile.Splitting/Snakefile  \
--jobs 1 \
--configfile $config_file \
--rerun-incomplete

wget https://github.com/mikefarah/yq/releases/download/v4.44.5/yq_linux_amd64.tar.gz -O - |\
tar xz 

chmod +x yq_linux_amd64 
output_jobs=$(cat $($(realpath yq_linux_amd64) '.mutect2_path' $config_file)/splits/JOBS)

echo $output_jobs


sbatch -J HEAD -p priopark -t 88:00:00 --mail-type=END --mem 3GB  --wrap "snakemake \
--snakefile /n/data1/hms/dbmi/park/dominika/testing/smaht/pipelines/MosaicForecast/Snakefile.MF/Snakefile  \
--executor slurm \
--default-resources slurm_account=park_contrib \
slurm_partition=park \
--jobs 1000 \
--rerun-incomplete \
--configfile \
$config_file \
--config jobs=$output_jobs"

rm yq_linux_amd64 
rm yq.1
