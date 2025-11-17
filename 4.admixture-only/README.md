# Admixture-Only Variant Identification and Visualization
The admixture only variants were broadly identified by intersecting the RUFUS and Mutect2 calls for the COLOBLT50 cell line with the truth set found in this publication. These data sets served as the inputs to all visualizations generated in Figures 3 and Figures S4.

## Creating admixture only vcfs with the union of RUFUS and Mutect2 calls
1. The RUFUS variant calling algorithm was called using low-frequency detection parameters in 1mb mode using a paired contrl. Each run utilized a Genome Characterization Center (GCC) -specific sequencing bam (e.g. BCM, Broad, NYGC, UW, WashU) as the subject file, and the COLO829 BL-only bam (i.e. `COLO829BL_Ill_230X.bam`) for the paired control file. The full script to launch all jobs can be found in `launch_split_runs.sh`, but an example region call looked like:
```bash
bash \$RUFUS_ROOT/runRufus.sh -s $SUBJECT -c $CONTROL -r $REFERENCE -f ${RESOURCE_DIR}GRCh38_full_analysis_set_plus    _decoy_hla.25.Jhash -m 5 -k 25 -t 40 -L -vs -R chr${curr_chr}:${start_coord}-${end_coord}
```

2. Each region vcf was then trimmed and combined using `trim_and_combine.sh`, and normalized using `post_process.sh`. These vcfs are found in the `rufus_vcfs` directory.
3. Mutect2 vcfs, called in a similar fashion of one GCC vs the COLO829 BL-only paired control, were procured from the Data Acquisition Center (DAC) and Peter Park's group at Harvard Medical School. These vcfs are found in the `mutect2_vcfs` directory.
4. Each Mutect2 and RUFUS vcfs specific to a GCC were intersected (e.g. `mutect2_vcfs/NORMED.PASS_SNPs_ONLY.COLO829BLT-BCM.vcf.gz` with `rufus_vcfs/sorted.atomed.normed.bcm.400x.merged.SMHTCOLO829BLT50.FINAL.no_inherited.SNV_only.hdaf_single_sample.vcf.gz`)
5. The unions of these intersection was concatenated into a single vcf.
6. The union vcfs were then intersected with bitmasks representing the possible degrees of intersection for the five GCC-sets (i.e. 1-5), and intersected with the DAC-published truth set of the COLOBLT50 mixture to retain the complement (i.e. those **not** identified in the truth set - `SMaHT_COLO829_SNV_truth_set_v1.0.vcf.gz`). This created five vcfs - one containing admixture-only variants that were found in a single GCC data set (i.e. BCM or Broad or NYGC or UW or WashU), one containing admixture-only variants found in two GCCs (i.e BCM and Broad or BCM and NYGC...), ... up to the vcf with admixture-only variants found in all five GCCs.
7. These vcfs from step 6 served as the input to the next visualization steps detailed below. The code describing steps 4-6 can be found in `mut_ruf_isec_all.sh`

## Figure 3 Generation

### Upset Plot
1. Mutect2 and RUFUS vcfs per GCC, plus the DAC-generated truth-set, were fed into `mut_ruf_isec_all.sh` to create all possible intersections of the call-sets that are not contained in the truth-set. 
2. The output from step 1 gets saved into a file `upset_counts.txt`. These counts were copied and pasted into the visualization script `upset_plot.ipynb` for Fig. 3C generation.

## AF Histogram Plot
1. Taking the vcfs from `mut_ruf_isec_all.sh`, allele frequencies were calculated using `batch_filter_afs.sh`, and then counted with `compose_counts.sh`.
2. Counts were manually transferred from the output of `composed_counts.sh` into arrays in `af_histo.ipynb`. The python script was then run to generate Fig. 3D.

## Supplemental Figure 4 Generation
1. Pileups were performed at each set of sites corresponding to the degree of intersection vcfs described in the top section of this README, using `pileup.slurm` for the subject BLT50 control bams, and `negative_ctrl_pileup.slurm` for the negative control bams.
2. Resulting pileup vcfs per subject/negative-control category were then merged into a single vcf. Single vcfs were then filtered for called-ALT entries only, and REF entries. Each ALT or REF vcf was then filtered by alternate observations or depth, respectively, and filtered by the number of observations. The final filtered vcf counts were tallied to a plain text file, counts.out. All of this was performed in `merge_and_tally.sh`
3. `transpose_counts.sh` was then run on the counts.out file to make copying the counts easily into the final python visualization script
4. `s4_double_histo.ipynb` and `s4_line_chart.ipynb` were then run with these counts to provide the histogram visualizations and line chart visualization seen in Figure S4.
