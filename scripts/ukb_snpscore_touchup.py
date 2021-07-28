#!/bin/python


import pandas as pd



raw_dir = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw'


# autism -- add columns
gwas_df = pd.read_table(os.path.join(raw_dir, 'daner_AUT_meta14_CEU_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tsv'), header=None)

gwas_df.columns = [
    'CHR',
    'BP',
    'SNP',
    'A1',
    'A2',
    'OR',
    'lb95',
    'ub95',
    'EFFECT',
    'SE',
    'P',
    'FRQ_A1',
    'INFO',
    'N',
    'DIRECTION'
]

gwas_df.to_csv(os.path.join(raw_dir, 'fmt_daner_AUT_meta14_CEU_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tsv'), index=None)


/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t%AF\t[%ES\t%SE\t%SS\t%NC]\n' \
./ieu-a-801.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tA1_FRQ\tBETA\tSE\tN\tNCases"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> BipolarDisorder_2011_NatGen_Sklar.txt


cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ieu-b-2.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> Alzheimers_2019_NatGen_Kunkle.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE\t%SS]\n' \
./ieu-a-812.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE\tN"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> Parkinsons_2009_NatGen_SimonSanchez.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE\t%AF]\n' \
./prot-a-2179.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE\tA1_FRQ"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> e3ParkinProtein_2018_Nat_Sun.txt


cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST005536.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> DiabetesType1_2015_NatGen_Onengut.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST001198.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> MultipleSclerosis_2011_Nat_Sawcer.txt


cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST005523.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> CeliacDisease_2011_NatGen_Trynka.txt


cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ieu-a-832.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> RheumatoidArthritus_2014_Nat_Okada.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST004131.vcf| \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> RheumatoidArthritus_2014_Nat_Okada.txt




cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE\t%AF]\n' \
./ebi-a-GCST003116.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE\tA1_FRQ"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> CoronaryArteryDisease_2015_NatGen_Nikpay.txt




cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE\t%AF]\n' \
./ebi-a-GCST000998.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE\tA1_FRQ"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> CoronaryArteryDisease_2011_NatGen_Schunkert.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST000755.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> ChronicKidneyDisease_2016_NatComm_Pattaro.txt




cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE\t%AF]\n' \
./ebi-a-GCST003045.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE\tA1_FRQ"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> UlcerativeColitis_2015_NatGen_Liu.txt




cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST001790.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> SerumGout_2013_NatGen_Kottgen.txt


cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST001791.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> SerumUrate_2013_NatGen_Kottgen.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST000755.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> HDLCholesterol_2010_Nat_Teslovich.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST000760.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> TotalCholesterol_2010_Nat_Teslovich.txt


cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST002920.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> Neuroticism_2015_JAMAPsych_deMoor.txt




cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ieu-a-29.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> BirthLength_2015_HumMolGen_vanderValk.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ieu-a-1096.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> ChildhoodObesity_2012_NatGen_Bradfield.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ieu-a-28.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> HeadCircumfrence_2012_NatGen_Taal.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST005647.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> AmyotrophicLateralSclerosis_2018_Neuron_Nicolas.txt


cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ieu-a-1108.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> IschemicStroke_2016_Neurology_Malik.txt



cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ieu-a-1108.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> CardioembolicStroke_2016_Neurology_Malik.txt


cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST006910.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> IschemicStroke_2018_NatGen_Malik.txt


cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw
/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/bin/bcftools \
query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t[%ES\t%SE]\n' \
./ebi-a-GCST006907.vcf | \
awk 'BEGIN {print "SNP\tP\tCHR\tPOS\tA1\tA2\tBETA\tSE"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' \
> LargeArteryAtherosclerosis_2018_NatGen_Malik.txt

