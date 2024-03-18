#this maps the fastq.gz files to reference fasta (rcrs.fa) this assumes that fastq.gz is single end sequenced or paired reads already merged and sequencing adapters removed (this is how sequencing data is usually, but not always published.
#If the input data is published whole genome data in the bam format I prefer to uses samtools bam2fq first and then remap the sequences to rcrs. It's beter then extracting the MT reads from whole genome data (this way you lose coverage in coservative regions).
bwa aln -l 16500 -n 0.01 -o 2 -t 4 rcrs.fa fastq.gz | bwa samse rcrs.fa - fastq.gz| samtools view -h -Su -F 4 - | samtools sort -T temp_file name -o output_file_name.bam

#removing sequence duplicates 2 approches
samtools rmdup -S output_file_name.bam output_file_name_unique.bam
#this requires python script from Mayer Kircher 2010 paper (https://github.com/shendurelab/cfDNA/blob/master/FilterUniqueBAM.py), and we use it in our oficial pipeline. However at some point samtools introduced thier own tool, and this seems to do the job, although I did not compared the two approaches statistically.
samtools view -F 4 -h output_file_name.bam | FilterUniqueSAMCons.py | samtools view -h -Su - > output_file_name_unique.bam

#This is optional and something we inherited from Uppsala. It removes reads shorter than 35bp and with lower than 90% identity with reference. This is a filter that only we (plus Uppsala, Stokholm and Ankara) are using and the python script is unpublished and I am not sure whatever I am allowed to share it so maybe omit this step.
samtools calmd output_file_name_unique.bam rcrs.fa | percidentity_threshold.py 0.9 35 short.txt | samtools view -bS - > output_file_name_unique_90prc.bam

#I am using two ways to build consensus. Using Samtools is something developed for a course we are teaching, as it does not require additional tools and it seems to do the job faster and with the same result as angsd. Stil angsd in our official pipeline and I did not have the time to systematically compare them.
samtools consensus -m simple -d 3 --min-MQ 30 --min-BQ 30 -c 0.5 output_file_name_unique.bam -o output_consensus_file_name.fa

angsd -doFasta 2 -doCounts 1 -i utput_file_name_unique.bam  -minQ 30 -minMapQ 30 -setMinDepth 3 -out output_consensus_file_name
#angsd by default produces gziped files
gunzip output_consensus_file_name.fa

#Either way you will have to replace the header of the fasta file so it would contain the sample name. On default it is ">MT" for samtools and ">" for angsd. I use something like
sed -i "s/>/>$(basename output_consensus_file_name.fa .fa)/" output_consensus_file_name.fa