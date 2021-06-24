#!/bin/bash
#example: psub.sh --page3 --study aric_ea --trait epicare_qtdur --file pheno_ecg 

sugen="/nas02/depts/epi/CVDGeneNas/sol/taor/SUGEN_github/SUGEN/SUGEN"

#store command for recursion
cmd="$0 $@"

#set relative program directories
pipe_dir="$(dirname "$0")"
this_path="$(realpath $0)"
pheno_dir=$pipe_dir/pheno
geno_dir=$pipe_dir/geno

#defaults
page="pageiii"
outdir="."

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
  *)
  echo "Unknown option $1"
  exit
	;;
esac

shift
done

#if trait is file containing traits, then run once for each line
if [[ -f $trait ]]; then
  this_script="$(realpath $0)"
  while read t; do
    [[ -z "$t" ]] && exit
    bash $cmd --trait $t
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
vcflist=$basedir/`yq eval .$study.vcflist $geno_dir/${page}_geno.yml`*.vcf.gz

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
echo "page: $page" >> $ymlout
echo "pheno: $pheno" >> $ymlout
echo "id: $id" >> $ymlout
echo "fid: $fid" >> $ymlout
echo "formula: $formula" >> $ymlout
echo "basedir: $basedir" >> $ymlout
echo "date: \"`date`\"" >> $ymlout

[[ "$test" == true ]] && exit

# cat $ymlout

echo "vcflist: " >> $ymlout

#loop over each vcf file
pidlist=""
for f in $vcflist; do
  fn=`basename $f`
  echo "$trait: Submitting $fn"
  sout=` sbatch -o $logsdir/log_${study}__${trait}.log --wrap="$sugen --pheno $pheno_file --id-col $id --family-col $fid --vcf $f --formula $formula --unweighted --out-prefix $resdir/${study}__${trait}__$fn --dosage" `
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