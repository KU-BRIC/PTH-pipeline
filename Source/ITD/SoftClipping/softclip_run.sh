#!/usr/bin/env bash



while getopts hi:o:s:g:m: opt; do
    case $opt in
    h)
        echo "Usage: $0 -i in_bam -o out_dir -s sample_name -g ref_gen -m min_cigar_spike [15]"
        exit
        ;;
    i)
        in_bam=$OPTARG
        echo "in_bam = $in_bam"
        ;;
    o)
        out_dir=$OPTARG
        echo "out_dir = $out_dir"
        ;;
    s)
	sample=$OPTARG
	echo "sample = $sample"
	;;
    g)
	ref_gen=$OPTARG
	echo "ref_gen = $ref_gen"
	;;
    m)
	min_cigar_spike=${OPTARG}
	echo "min_cigar_spike = $min_cigar_spike"
	;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        ;;
    esac
done


# set default value if not set
min_cigar_spike=${min_cigar_spike:-15}
echo "min_cigar_spike = $min_cigar_spike"


function soft_clip_spikes () {
        local in_bam="$1"
        local softclip_dir="$2"
        local sample="$3"
	local REF_GEN=$ref_gen

	if [ ! -d $softclip_dir ]; then 
		mkdir -p $softclip_dir
	fi

        local sample_chr_pos_cigar_seq=$softclip_dir/${sample}_spcs.tsv

        # extract pos cigar seq from BAM from the FLT3 region
        echo "Extracting FLT3 data from $sample"
        FLT3ITD="chr13:28033800-28034200"
        #samtools view $in_bam $FLT3ITD | cut -f3,4,6,10 | awk -v s=$sample -v OFS="\t" -v refg=$REF_GEN '{print s, $1, $2, $3, $4}' > $sample_chr_pos_cigar_seq
samtools view $in_bam $FLT3ITD | cut -f3,4,6,10 | awk -v s=$sample -v OFS="\t" '{print s, $1, $2, $3, $4}' | while read s r p c q; do echo -ne "$s\t$r\t$p\t"; samtools faidx $REF_GEN $r:$p-$p | tail -1 | tr -d "\n"; echo -e "\t.\t$c\t$q"; done > $sample_chr_pos_cigar_seq

        # call spike finding R script
        echo "Finding SoftClip Spikes for $sample"
        MIN_CIGAR_SPIKE=15
        Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/SoftClipping/softclip_spikes.R $sample $softclip_dir $sample_chr_pos_cigar_seq $MIN_CIGAR_SPIKE $REF_GEN
}


soft_clip_spikes $in_bam $out_dir $sample $ref_gen
