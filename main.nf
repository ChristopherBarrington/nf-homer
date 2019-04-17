


// conda environment to include
// Homer
// Juicer tools
// R
// 	tidyverse
// bowtie2

// paths that will be dynamic
BOWTIE2_INDEX='/camp/svc/reference/Genomics/iGenomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome'
GENOME_FASTA='/camp/svc/reference/Genomics/iGenomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa'
GENOME_FAI='/camp/svc/reference/Genomics/iGenomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa.fai'

// get the input FastQ files into a <SAMPLE><LIBRARY><LANE><FASTQ> format
Channel
	.from(params.raw_data.samples.collect{
		sample,
		libraries -> libraries.collect{
			library,
			lanes -> lanes.collect{
				lane,
				files -> [sample,library,lane,files]}}.sum()}.sum())
	.transpose()
	.map{ v -> [v[0], v[1], v[2], (v[3] =~ '.*_(R.).*')[0][1], file(v[3])] }
	.dump(tag: 'INPUT_fastq_files_to_subset')
	.set{INPUT__fastq_files_to_subset}


/***********************

 Preprocess FastQ files

 ***********************/


// if the --split parameter is used, split the FastQ files into bitesize chunks before continuing
process split_input_files_into_chunks {

	tag {[sample, library, lane, read_pair].join('-')}

	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'

	input:
		set val(sample), val(library), val(lane), val(read_pair), file(fastq) from INPUT__fastq_files_to_subset

	output:
		set val(sample), val(library), val(lane), val(read_pair), file('*fastq') into OUTPUT__split_input_files_into_chunks

	script:
		if(params.raw_data.get('split_input_files', false)) {
			n_records = params.raw_data.get('split_input_files').toString().isInteger() ? params.raw_data.get('split_input_files') : 500000
			n_lines = n_records * 4
			"""
			gunzip --to-stdout ${fastq} | split --lines $n_lines --suffix-length 5 --additional-suffix .fastq - ''
			"""
		} else {
			"""
			gunzip --force --to-stdout ${fastq} > aaaaa.fastq
			"""
		}
}

OUTPUT__split_input_files_into_chunks
	.transpose()
	.map{ v -> [v[0], v[1], v[2], v[3], (v[4].toString() =~ '^.*/(.*?).fastq$')[0][1], file(v[4])] }
	.dump(tag: 'OUTPUT__split_input_files_into_chunks')
	.set{INPUT__trim_restriction_sites_from_fastq}

// read each input FastQ file and remove the restriction site (if present)
process trim_restriction_sites_from_fastq {

	tag {[sample, library, lane, read_pair, chunk].join('-')}

	module 'Homer/4.10-Perl-5.26.1-foss-2018a'

	input:
		set val(sample), val(library), val(lane), val(read_pair), val(chunk), file(fastq) from INPUT__trim_restriction_sites_from_fastq

	output:
		set val(sample), val(library), val(lane), val(read_pair), val(chunk), file('*.fastq.trimmed') into OUTPUT__trim_restriction_sites_from_fastq__TRIMMED_FASTQ

	script:

		"""
		homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 $fastq
		"""
}

OUTPUT__trim_restriction_sites_from_fastq__TRIMMED_FASTQ
	.dump(tag: 'OUTPUT__trim_restriction_sites_from_fastq__TRIMMED_FASTQ')
	.set{INPUT__align_trimmed_fastq}


/*****************************************************

 Align trimmed read pairs independently to the genome

 *****************************************************/


process align_trimmed_fastq {

	tag {[sample, library, lane, read_pair].join('-')}
	publishDir mode: 'copy', overwrite: true, path: "output/bowtie2/logs", pattern: '*log'

	module 'Bowtie2/2.3.4.1-foss-2018a'
	cpus 1
	time '00:30:00'
	memory '1G'
	executor 'slurm'

	input:
		set val(sample), val(library), val(lane), val(read_pair), val(chunk), file(fastq) from INPUT__align_trimmed_fastq

	output:
		set val(sample), val(library), val(lane), val(read_pair), val(chunk), file('*.sam') into OUTPUT__align_trimmed_fastq__ALIGNED_SAM
		file(output_log)

	script:
		output_file_root = [sample, library, lane, fastq.toString().replaceAll('.fastq.trimmed$', ''), read_pair].join('-')
		output_sam = output_file_root + '.sam'
		output_log = output_file_root + '.log'
		"""
		bowtie2 -p ${task.cpus} -x $BOWTIE2_INDEX -U $fastq > $output_sam
		ln -s .command.err $output_log
		"""
}

OUTPUT__align_trimmed_fastq__ALIGNED_SAM
	.dump(tag: 'OUTPUT__align_trimmed_fastq__ALIGNED_SAM')
	.into{INPUT__convert_sam_chunks_to_bam;
	      INPUT__make_tag_directory}

// convert the mapped sam chunks into bam
////process convert_sam_chunks_to_bam {
////	
////	module 'SAMtools/1.9-foss-2018b'
////
////	input:
////		set val(sample), val(library), val(lane), val(read_pair), val(chunk), file(sam) from INPUT__convert_sam_chunks_to_bam
////
////	output:
////		set val(sample), val(library), val(lane), val(read_pair), val(chunk), file('*.bam') into OUTPUT__convert_sam_chunks_to_bam__BAM
////
////	script:
////		output_file = sam.toString().replaceAll('.sam$', '.bam')
////		"""
////		samtools view --output-fmt=BAM -o $output_file $sam
////		"""
////}
////
////OUTPUT__convert_sam_chunks_to_bam__BAM
////	.dump(tag: 'OUTPUT__convert_sam_chunks_to_bam__BAM')
////	.set{INPUT__merge_map_chunks}

// combine the mapped fragments from lane chunks
////process merge_map_chunks {
////
////	module 'SAMtools/1.9-foss-2018b'
////
////	input:
////		set val(sample), val(library), val(lane), val(read_pair), val(chunk), file(sams) from INPUT__merge_map_chunks.groupTuple(by: [0,1,2,4]).dump(tag: 'INPUT__merge_map_chunks')
////
////	output:
////		set val(sample), val(library), val(lane), val(read_pair), val(chunk), file(output_file) into OUTPUT__merge_map_chunks_BAM
////
////	script:
////		output_file = read_pair + '.sam'
////		"""
////		SAMS=`echo $sams | tr ' ' '\n' | sort`
////		samtools cat \${SAMS[@]} | samtools view -h -o $output_file -
////		"""
////}
////
////OUTPUT__merge_map_chunks_BAM
////	.dump(tag: 'OUTPUT__merge_map_chunks_BAM')


/******************************************************

 Make tag directories for lanes, libraries and samples
  
 ******************************************************/


// make a tag directory, grouping the chunks together into lane-level datasets
// PCR duplicates removed using the `-tbl 1` argument
INPUT__make_tag_directory
	.groupTuple(by: [0,1,2,4])
	.groupTuple(by: [0,1,2])
	.map{ v -> [v[0], v[1], v[2], v[3].flatten().sort(), v[4].sort(), v[5].flatten().sort{ it.getFileName() }]}
	.dump(tag: 'INPUT__make_tag_directory_for_lanes__LANE_GROUPS')
	.set{INPUT__make_tag_directory_for_lanes__LANE_GROUPS}

process make_tag_directory_for_lanes {

	tag {[sample, library, lane].join('-')}
	publishDir mode: 'copy', overwrite: true, path: 'output/tags', pattern: "$output_dir/*tsv"

	module 'Homer/4.10-Perl-5.26.1-foss-2018a'
	cache true
	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'
	
	input:
		set val(sample), val(library), val(lane), val(read_pairs), val(chunks), file(sams) from INPUT__make_tag_directory_for_lanes__LANE_GROUPS
		// `read_pairs` and `chunks` not used; could remove v[3] and v[4]

	output:
		set val(sample), val(library), val(lane), file(output_dir) into OUTPUT__make_tag_directory_for_lanes
		set val(output_dir), file(output_dir) into OUTPUT__make_tag_directory_for_lanes__QC, OUTPUT__make_tag_directory_for_lanes__CONVERT_FORMAT

	script:
		output_dir = [sample, library, lane].join('-')
		"""
		SAMS=`echo $sams | tr ' ' '\n' | paste --delimiters=,  - -`
		makeTagDirectory $output_dir -genome $GENOME_FASTA -restrictionSite GATC -removePEbg -tbp 1 \${SAMS[@]}
		"""
}

OUTPUT__make_tag_directory_for_lanes
	.groupTuple(by: [0,1])
	.map{ v -> [v[0], v[1], v[2].flatten().sort(), v[3].flatten().sort{ it.getFileName() }]}
	.dump(tag: 'INPUT__merge_lane_tags_into_libraries')
	.set{INPUT__merge_lane_tags_into_libraries}


// merge the lane-level datasets into library-level datasets
// PCR duplicates removed using the `-tbl 1` argument
process merge_lane_tags_into_libraries {

	tag {[sample, library].join('-')}
	publishDir mode: 'copy', overwrite: true, path: 'output/tags', pattern: '**tsv'
	publishDir mode: 'copy', overwrite: true, path: 'output/stats', pattern: '**txt'

	module 'Homer/4.10-Perl-5.26.1-foss-2018a'
	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'

	when:
		true

	input:
		set val(sample), val(library), val(lanes), file(lane_tag_dirs) from INPUT__merge_lane_tags_into_libraries
		// could drop `lanes` here amd count `lane_tag_dirs` and remove v[2]

	output:
		set val(sample), val(library), file(output_dir) into OUTPUT__merge_lane_tags_into_libraries
		set val(output_dir), file(output_dir) into OUTPUT__merge_lane_tags_into_libraries__QC, OUTPUT__merge_lane_tags_into_libraries__CONVERT_FORMAT, OUTPUT__merge_lane_tags_into_libraries__ANALYZE_HIC
		file "$output_dir/*{tsv,txt}" into OUTPUT__merge_lane_tags_into_libraries__TO_PUBLISH mode flatten

	script:
		output_dir = [sample, library].join('-')
		if( [lanes].flatten().size() > 1 ) {
			"""
			makeTagDirectory $output_dir -genome $GENOME_FASTA -restrictionSite GATC -tbp 1 -d $lane_tag_dirs
			"""
		} else {
			"""
			ln -s $lane_tag_dirs $output_dir
			"""
		}
}

OUTPUT__merge_lane_tags_into_libraries__TO_PUBLISH
	.dump(tag: 'OUTPUT__merge_lane_tags_into_libraries__TO_PUBLISH')

OUTPUT__merge_lane_tags_into_libraries
	.groupTuple(by: 0,)
	.map{ v -> [v[0], v[1].flatten().sort(), v[2].flatten().sort{ it.getFileName() }]}
	.dump(tag: 'OUTPUT__merge_lane_tags_into_libraries')
	.set{INPUT__merge_library_tags_into_samples}

// merge the lane-level datasets into library-level datasets
// PCR duplicates are not removed here, since these are biological replicates
process merge_library_tags_into_samples {

	tag {sample}
	publishDir mode: 'copy', overwrite: true, path: 'output/tags', pattern: '**tsv'
	publishDir mode: 'copy', overwrite: true, path: 'output/stats', pattern: '**txt'

	module 'Homer/4.10-Perl-5.26.1-foss-2018a'
	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'

	when:
		true

	input:
		set val(sample), val(libraries), file(library_tag_dirs) from INPUT__merge_library_tags_into_samples
		// could drop `libraries` here amd count `library_tag_dirs` and remove v[1]

	output:
		set val(sample), file(sample) into OUTPUT__merge_library_tags_into_samples
		set val(sample), file(sample) into OUTPUT__merge_library_tags_into_samples__QC, OUTPUT__merge_library_tags_into_samples__CONVERT_FORMAT, OUTPUT__merge_library_tags_into_samples__ANALYZE_HIC
		file "$sample/*{tsv,txt}" into OUTPUT__merge_library_tags_into_samples__TO_PUBLISH mode flatten

	script:
		if( [libraries].flatten().size() > 1 ) {
			"""
			makeTagDirectory $sample -genome $GENOME_FASTA -restrictionSite GATC -d $library_tag_dirs
			"""
		} else {
			"""
			ln -s $library_tag_dirs $sample
			"""
		}
}

OUTPUT__merge_library_tags_into_samples
	.dump(tag: 'OUTPUT__merge_library_tags_into_samples')

OUTPUT__merge_library_tags_into_samples__TO_PUBLISH
	.dump(tag: 'OUTPUT__merge_library_tags_into_samples_TO_PUBLISH')


/**********************************************
 
 Convert tag directories to other Hi-C formats
 
 **********************************************/

OUTPUT__make_tag_directory_for_lanes__CONVERT_FORMAT
	.concat(OUTPUT__merge_lane_tags_into_libraries__CONVERT_FORMAT, OUTPUT__merge_library_tags_into_samples__CONVERT_FORMAT)
	.dump(tag: 'INPUT__convet_format_processes')
	.into{INPUT__make_juicebox_file_from_maps;
	      INPUT__make_distiller_file_from_maps}

process make_juicebox_file_from_maps {

	tag {dataset_name}
	publishDir mode: 'copy', overwrite: true, path: 'output/Juicebox'

	module 'Homer/4.10-Perl-5.26.1-foss-2018a'
	cache true
	cpus 1
	time '00:30:00'
	memory '1G'
	executor 'local'

	when:
		false

	input:
		set val(dataset_name), file(tag_directory_path) from INPUT__make_juicebox_file_from_maps

	output:
		set val(dataset_name), file(output_file) into OUTPUT__make_juicebox_file_from_maps

	script:
		output_file = dataset_name + '.hic'
		resolutions = params.hic_file_options.get('resolutions', [2000,5000,10000,25000,50000,100000,200000]).sort().join(',')
		"""
		# convert the tag.tsv files into 'Juicer short format' using dummy fragment numbers (0/1)
		# convert the sorted short-format file into .hic
		tagDir2hicFile.pl $tag_directory_path/ -juicer $output_file -genome <(cut -f 1,2 $GENOME_FAI) -juicerExe /camp/stp/babs/working/barrinc/Git-clones/juicer/SLURM/scripts/juicer_tools -juicerOpt '-r $resolutions' -p ${task.cpus}
		"""
}


/***************************************

Make QC plots for libraries and samples

 ***************************************/


// plot QC metrics for the datasets
OUTPUT__make_tag_directory_for_lanes__QC
	.concat(OUTPUT__merge_lane_tags_into_libraries__QC, OUTPUT__merge_library_tags_into_samples__QC)
	.dump(tag: 'INPUT__plot_tag_qc')
	.into{INPUT__plot_tag_qc_petagLocalDistribution;
	      INPUT__plot_tag_qc_petagFreqDistribution;
	      INPUT__plot_tag_qc_petagRestrictionDistribution}

process plot_tag_qc_petagLocalDistribution {

	tag {dataset_name}
	publishDir mode: 'copy', overwrite: true, path: 'output/qc/petagLocalDistribution'

	module 'UDUNITS/2.2.26-foss-2016b'
	module 'R/3.5.1-foss-2016b-BABS'
	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'
	cache true

	when:
		true

	input:
		set val(dataset_name), file(tag_directory_path) from INPUT__plot_tag_qc_petagLocalDistribution

	output:
		set val(dataset_name), file("*$output_format") into OUTPUT__plot_tag_qc_petagLocalDistribution

	script:
		output_format = 'pdf'
		template 'plot_petagLocalDistribution.R'
}

process plot_tag_qc_petagFreqDistribution {

	module 'UDUNITS/2.2.26-foss-2016b'
	module 'R/3.5.1-foss-2016b-BABS'
	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'
	cache true

	tag {dataset_name}
	publishDir mode: 'copy', overwrite: true, path: 'output/qc/petagFreqDistribution'

	when:
		true

	input:
		set val(dataset_name), file(tag_directory_path) from INPUT__plot_tag_qc_petagFreqDistribution

	output:
		set val(dataset_name), file("*$output_format") into OUTPUT__plot_tag_qc_petagFreqDistribution

	script:
		output_format = 'pdf'
		template 'plot_petagFreqDistribution.R'
}

process plot_tag_qc_petagRestrictionDistribution {

	module 'UDUNITS/2.2.26-foss-2016b'
	module 'R/3.5.1-foss-2016b-BABS'
	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'
	cache true

	tag {dataset_name}
	publishDir mode: 'copy', overwrite: true, path: 'output/qc/petagRestrictionDistribution'

	when:
		true

	input:
		set val(dataset_name), file(tag_directory_path) from INPUT__plot_tag_qc_petagRestrictionDistribution

	output:
		set val(dataset_name), file("*$output_format") into OUTPUT__plot_tag_qc_petagRestrictionDistribution

	script:
		output_format = 'pdf'
		restriction_site = 'GATC'
		mismatches = 0
		template 'plot_petagRestrictionDistribution.R'
}


/****************************

 Run Homer tools on the data

 ****************************/


OUTPUT__merge_lane_tags_into_libraries__ANALYZE_HIC
	.concat(OUTPUT__merge_library_tags_into_samples__ANALYZE_HIC)
	.dump(tag: 'INPUT__analyze_hic')
	.into{INPUT__analyze_hic_create_matrix;
	      INPUT__analyze_hic_chromatin_compaction;
	      INPUT__analyze_hic_chromatin_compartments; // add PCA analysis from http://homer.ucsd.edu/homer/interactions2/HiCpca.html
	      INPUT__analyze_hic_find_tads_and_loops}

// make an interaction matrix from the tag directory
process analyze_hic_create_matrix {

	tag {dataset_name}
	publishDir mode: 'copy', overwrite: true, path: "output/interaction_matrices/$normalisation_method/$resolution"

	module 'Homer/4.10-Perl-5.26.1-foss-2018a'
	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'
	cache true

	resolutions = [2000, 10000, 50000, 100000] // get from params
	normalisation_methods = ['coverageNorm','distNorm']

	when:
		true

	input:
		set val(dataset_name), file(tag_directory_path) from INPUT__analyze_hic_create_matrix
		each resolution from resolutions
		each normalisation_method from normalisation_methods

	output:
		set val(dataset_name), val(normalisation_method), val(resolution), val(window), file(output_file) into OUTPUT__analyze_hic_create_matrix

	script:
		window = resolution
		output_file = [dataset_name, 'homer'].join('.')
		"""
		analyzeHiC $tag_directory_path -cpu ${task.cpus} -res $resolution -window $window -$normalisation_method -balance -o $output_file
		"""
}

// calculate distal-to-local (DLR) and interchromosomal fraction (ICF)
process analyze_hic_chromatin_compaction {

	tag {dataset_name}
	publishDir mode: 'copy', overwrite: true, path: "output/chromatin_compaction/$normalisation_method/res_$resolution/win_$window", pattern: '*bedGraph'
	publishDir mode: 'copy', overwrite: true, path: "output/chromatin_compaction/$normalisation_method/res_$resolution/win_$window/logs", pattern: '*log'

	module 'Homer/4.10-Perl-5.26.1-foss-2018a'
	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'
	cache true

	resolutions = [5000, 10000] // get from params
	normalisation_methods = ['coverageNorm']

	when:
		!params.containsKey('chromatin_compaction')

	input:
		set val(dataset_name), file(tag_directory_path) from INPUT__analyze_hic_chromatin_compaction
		each resolution from resolutions
		each normalisation_method from normalisation_methods

	output:
		set val(dataset_name), val(normalisation_method), val(resolution), val(window), file("$output_file*") into OUTPUT__analyze_hic_chromatin_compaction

	script:
		max_distance = 3000000
		window = resolution*2
		output_file = dataset_name
		output_log = dataset_name + '.log'
		"""
		analyzeHiC $tag_directory_path \
			-cpu ${task.cpus} \
			-res $resolution \
			-window $window \
			-$normalisation_method \
			-nomatrix \
			-compactionStats $output_file \
			-dlrDistance $max_distance
		ln -s .command.err $output_log
		"""
}

process analyze_hic_find_tads_and_loops {

	tag {dataset_name}
	publishDir mode: 'copy', overwrite: true, path: "output/loops_and_tads/res_$resolution/win_$window", pattern: '*{bedGraph,txt,bed}'
	publishDir mode: 'copy', overwrite: true, path: "output/loops_and_tads/res_$resolution/win_$window/logs", pattern: '*log'

	module 'Homer/4.10-Perl-5.26.1-foss-2018a'
	cpus 1
	time '00:30:00'
	memory '6G'
	executor 'slurm'
	cache true

	resolutions = [5000, 10000] // get from params

	when:
		!params.containsKey('chromatin_compaction')

	input:
		set val(dataset_name), file(tag_directory_path) from INPUT__analyze_hic_find_tads_and_loops
		each resolution from resolutions

	output:
		set val(dataset_name), val('none'), val(resolution), val(window), file("$output_file*") into OUTPUT__analyze_hic_find_tads_and_loops

	script:
		window = resolution*2
		min_tad_size = 5000
		min_tad_score = 1.5
		keep_overlapping_tads = true
		min_distance = 2000
		max_distance = 2000000
		insulation_distance = 10000
		directionality_index_distance = 30000
		output_file = dataset_name
		output_log = dataset_name + '.log'
		"""
		# can include -p flag here to blacklist regions (launches the filterTADsAndLoops.pl script)
		findTADsAndLoops.pl find $tag_directory_path \
			-cpu ${task.cpus} \
			-res $resolution \
			-window $window \
			-balance \
			-minDist $min_distance \
			-maxDist $max_distance \
			-insDist $insulation_distance \
			-diDist $directionality_index_distance \
			-o $dataset_name \
			-override
		ln -s .command.err $output_log
		"""
}
