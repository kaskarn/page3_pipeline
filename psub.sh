#!/bin/bash
#example: psub.sh --page3 --study aric_ea --trait epicare_qtdur --file pheno_ecg 

sugen="/nas02/depts/epi/CVDGeneNas/sol/taor/SUGEN_github/SUGEN/SUGEN"

#store command for recursion
cmd="$0 $@"

#set relative program directories

this_path="$(realpath $0)"
this_dir="$(dirname $this_path)"
pheno_dir=$this_dir/pheno
geno_dir=$this_dir/geno

#defaults
page="pageiii"
outdir="."

#halp
if [ $# -lt 1 ]
then
  echo -e '

  SUGEN / SLURM wrapper for PAGE GWAS
  (Antoine Baldassari -- baldassa@email.unc.edu)
  
  USAGE
  bash psub.sh --study studyname --pheno phenoname --trait traitname [--outdir out --covars varlist --test --pageii]

  REQUIRED ARGUMENTS
  --study STUDY : Study name  
  One of: whims, garnet, gecco, hipfx, mopmap, aric_aa, aric_ea, mega_all, mega_aa, mega_ha

  --trait TRAIT : Single trait name, or file listing one trait name per line
  Example: "--trait qt_duration"; "--trait mytraits.txt"

  --pheno PHEN : Phenotype file, or YAML specification file (with, or without the .yml extension)
  Phenotype file specifications are yml-formatted files contained in ./pheno,
  relative to this script. These files specify covariate adjustments,
  id names, study- or trait- specific covariates, trait models, and others. See ./page_ecg.yml
  for an example.
  Examples: 
    --pheno page_ecg ... (use Antoine'"'"'s ECG phenotype file, described in ./page_ecg.yml)

  OPTIONAL ARGUMENTS
  -h | --help:  print this manual and exit

  --covars COVARS: quoted list of space-delimited covariates overriding phenotype specification
  Example: --covars "age bmi pc1 pc2"

  --pageii : Use PAGE II genotypes instead of PAGE III
  --test : Set up only; does not run SUGEN
  --test1 : Only run the first CHR

'
  exit
fi

#parsing
while [[ $# -gt 0 ]]
do
case $1 in
  --pageii)
  page="pageii"
  ;;
  --pageiii)
  page="pageiii"
  ;;
	-s|--study)
	study="$2"
	shift
	;;
	-t|--trait)
	trait="$2"
	shift
	;;
	-c|--covars)
	covars_cmd="$2"
	shift
	;;
  -p|--pheno)
  pheno="$2"
  shift
  ;;
  -o|--outdir)
  outdir="$2"
  shift
  ;;
  --test)
  test=true
  outdir="test"
  ;;
  --test1)
  test1=true
  outdir="test"
  ;;
  -h|--help)
  bash $0
  exit
  ;;
  *)
  echo "Unknown option $1"
  exit
	;;
esac

shift
done

#debugging output
# if [[ ! -z "$test" ]]; then
#   echo $this_path
#   echo $this_dir
#   echo $geno_dir
#   echo $pheno_dir
# fi

#if trait is file containing traits, then run once for each trait within
if [[ -f $trait ]]; then
  this_script="$(realpath $0)"
  while read t; do
    [[ -z "$t" ]] || bash $cmd --trait $t
  done < $trait
  exit
fi

# Testing
# study="aric_ea"
# page="pageiii"
# trait="epicare_qtdur"
# pheno="page_ecg"

#Set phenotype or YML file
if [[ -f $pheno ]]; then
  pheno_file=$pheno
  #error checking
else
  pheno=$pheno_dir/$pheno.yml
  pheno_file=` yq eval .path $pheno`
fi

#set covariates
if [[ -z $covars_cmd ]]; then

  covars_s=` yq eval .studies.$study.covars[] $pheno `
  covars_f=` yq eval .covars[] $pheno `
  covars_t=` yq eval .traits.$trait.covars[] $pheno`
  pc=` yq eval .pc_prefix $pheno`
  covars_pc=`for i in {1..10}; do echo $pc$i; done`
  covars=""
  [[ "$covars_s" == "null" ]] || covars="$covars $covars_s"
  [[ "$covars_f" == "null" ]] || covars="$covars $covars_f"
  [[ "$covars_t" == "null" ]] || covars="$covars $covars_t"
  [[ "$pc" == "null" ]] || covars="$covars $covars_pc"

else
  covars=$covars_cmd
fi

#set model
model="linear"
model_y=` yq eval .traits.$trait.model $pheno`
[[ "$model_y" == "null" ]] || model=$model_y


#set participant id column
if [[ -z $id_cmd ]]; then
  id=` yq eval .$page.id $pheno`
else
  id=$id_cmd
fi

#set family id column
if [[ -z $fid_cmd ]]; then
  fid=` yq eval .$page.fid $pheno`
else
  fid=$fid_cmd
fi

#set formula
eqn_rhs=`echo $covars | sed 's/ \+/+/g'`
formula="${trait}=$eqn_rhs"

#set list of vcf files to loop through
basedir=`yq eval .path $geno_dir/${page}_geno.yml `
vcflist=$basedir/`yq eval .studies.$study.vcflist $geno_dir/${page}_geno.yml`*.vcf.gz

#set output location
dir_out="$outdir/$study/$trait"
logsdir="$dir_out/slurmlogs"
resdir="$dir_out/res"
[[ -d $outdir ]] || mkdir $outdir
[[ -d $outdir/$study ]] || mkdir $outdir/$study
[[ -d $dir_out ]] || mkdir $dir_out
[[ -d $logsdir ]] || mkdir $logsdir
[[ -d $resdir ]] || mkdir $resdir


#Start yaml summary
ymlout="$dir_out/summary.yml"
echo "#This file describes analyses submitted from this directory, and can be used as pipeline input to recreate them" > $ymlout
echo "#See pipeline documentation" >> $ymlout
echo "cmd: "'"'$cmd'"' >> $ymlout
echo "page: $page" >> $ymlout
echo "pheno: $pheno" >> $ymlout
echo "id: $id" >> $ymlout
echo "fid: $fid" >> $ymlout
echo "formula: $formula" >> $ymlout
echo "basedir: $basedir" >> $ymlout
echo "date: \"`date`\"" >> $ymlout
echo "model: $model" >> $ymlout

[[ "$test" == true ]] && exit

# cat $ymlout

echo "vcflist: " >> $ymlout

#loop over each vcf file
pidlist=""
for f in $vcflist; do
  fn=`basename $f`
  echo "$trait: Submitting $fn"
  sout=` sbatch -o $logsdir/log_${study}__${trait}.log -t4-0 --mem=4GB --wrap="$sugen --pheno $pheno_file --id-col $id --family-col $fid --vcf $f --formula $formula --unweighted --out-prefix $resdir/${study}__${trait}__$fn --dosage" `
  pid=`grep -o '[0-9]\+$' <(echo $sout)`
  pidlist="$pidlist $pid"
  echo "  - $f" >> $ymlout

  [[ "$test1" == true ]] && break
done

#add pids
echo "pidlist: " >> $ymlout
for pid in $pidlist; do
  echo "  - $pid" >> $ymlout
done