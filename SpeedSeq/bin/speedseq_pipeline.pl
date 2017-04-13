#!/usr/bin/env perl

### speedseq_pipeline.pl ##############################################################################
# Run the SpeedSeq pipeline.

### INCLUDES ######################################################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Path::Class;
use File::Path;
use File::Spec;
use File::Basename;
use YAML qw(LoadFile);

#use NGS::Tools::SpeedSeq;
use lib '/projects/ccrstaff/jw24/software/SpeedSeq/lib/NGS/Tools';

use SpeedSeq;

### COMMAND LINE DEFAULT ARGUMENTS ################################################################
our %opts = (
    debug                      => 'N',
    speedseq_version               => '0.0.1',  # this software is not versioned, so just make one for now
    align_seq                  => 'Y',
    merge_bams                 => 'Y',
    call_simple_variants           => 'Y',
    cal_structural_variants       => 'Y',
    ref                        => '/util/ccr/data/Broad/gatk_bundle/2.8/human_g1k_v37_decoy.fasta',
	exclude_intervals          => '/util/common/bioinformatics/speedseq/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed',
    modules                    => 'samtools,root,speedseq,perl/5.20.2,python/anaconda'
);

GetOptions(
    \%opts,
    "help|?",
    "man",
    "speedseq_version|v:s" => \$opts{'speedseq_version'},
    "run_id:s",
    "config|c:s",
    "output_dir:s",
    "data_type:s",
    "debug:s"                      => \$opts{'debug'},
    "ref|r:s"                      => \$opts{'ref'},
    "modules:s"                    => \$opts{'modules'},
    "align_seq:s" => \$opts{'align_seq'},
    "call_simple_variants:s"    => \$opts{'call_simple_variants'},
    "call_structural_variants"             => \$opts{'call_structural_variants'},
    "merge_bams:s"       => \$opts{'merge_bams'}
) or pod2usage(2);

GetOptions();

####################################################################################################
# Check commandline sanity
####################################################################################################
if ( !$opts{'config'} or !$opts{'output_dir'} or !$opts{'run_id'} or !$opts{'data_type'}) {
    pod2usage(1);
}

if ( $opts{'help'} ) { pod2usage(1) }
if ( $opts{'man'} ) { pod2usage( -exitstatus => 0, -verbose => 2 ) }

my @list_of_modules = split( /\,/, $opts{'modules'} );

# Load the config yaml file
my $config;
if ( $opts{'config'} ) {
    $config = LoadFile( $opts{'config'} );
}
else {
    croak "Config yaml file not found. Please check your commandline\n";
}

####################################################################################################
# Root directories to store output files
####################################################################################################
my $parent_dir    = File::Spec->rel2abs( $opts{'output_dir'} );

my $speedseq = NGS::Tools::SpeedSeq->new(
    config        => $config,
	parent_dir    => $parent_dir,
    opts          => \%opts,
    modules       => \@list_of_modules
);

####################################################################################################
# Prepare directories for output files
####################################################################################################
# Switch to working directory
chdir($parent_dir);

# Initialize directory structure
$speedseq->dir_structure();

foreach my $phasedir ( 'ALN', 'SNV', 'SV' ) {
     my $phaseroot = $parent_dir . '/' . $phasedir;
     File::Path::make_path( File::Spec->rel2abs($phaseroot) );
     chdir( File::Spec->rel2abs($phaseroot) ) or croak "$!";
 
     my @phase_dirs = values %{ $speedseq->dir->{$phasedir} };
     File::Path::make_path(@phase_dirs);
 
     # Switch back to parent dir
     chdir( File::Spec->rel2abs($parent_dir) ) or croak "$!";
}

####################################################################################################
# Align raw sequences to ref genome
####################################################################################################

my $align_seq_ = $speedseq->submit_jobs('align');

my $merge_bams_ = $align_seq_->submit_jobs('merge_bam');

my $sv_ = $merge_bams_->submit_jobs('sv');


####################################################################################################
# Simple variant calling
####################################################################################################
# if we reach here, everything is ok
print "All jobs submitted...\n";

exit;

__END__


=head1 NAME

speedseq_pipeline.pl

=head1 SYNOPSIS

B<speedseq_pipeline.pl> <--config> <--run_id> <--output_dir> <--data_type>

	Arguments:
	--config			A YAML file containing sample names, lane names and paired fastq files for processing
	--run_id			A name/ID for this speedseq run, suggest to use the project, dataset or cohort name
	--output_dir			Output directory (for example: /scratch2/user/your_user_name).
        --data_type             	Indicate if the input sequence data are from whole genome sequencing or whole exome sequencying. valid values are "wgs" or "wes"
	
	Options:
	--help					Brief help message
	--man					Full documentation

	Examples:
	
	### cat config_speedseq.yaml
	#	SCVD0652:
	#	  CGCTCATT_H5LNFCCXX_L001:
	#	      R1: /projects/academic/big/BIG-Pilot-Projects/102515-MS-Study/Project_NOW_10587_B01_GRM_WGS.2015-09-19/Sample_CVD0652/fastq/CVD0652_CGCTCATT_H5LNFCCXX_L001_001.R1.fastq.gz
	#	      R2: /projects/academic/big/BIG-Pilot-Projects/102515-MS-Study/Project_NOW_10587_B01_GRM_WGS.2015-09-19/Sample_CVD0652/fastq/CVD0652_CGCTCATT_H5LNFCCXX_L001_001.R2.fastq.gz

	### Simple case, assuming everying will run smoothly and the default values are fine to you:
	perl /projects/ccrstaff/jw24/software/SpeedSeq/bin/speedseq_pipeline.pl --config config_speedseq.yaml --run_id MS --output_dir /gpfs/scratch/jw24/ --data_type wgs
	
	
=head1 AUTHOR

Jianxin Wang

Center for Computational Research, University at Buffalo

=head1 SEE ALSO

=cut

