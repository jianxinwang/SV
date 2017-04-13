package NGS::Tools::SpeedSeq;

use strict;
use warnings;
use Moose;
use 5.006;
use Carp;
use File::Basename;
use Params::Validate qw(:all);
use CCR::Utilities::General;

=head1 NAME

NGS::Tools::SpeedSeq - SpeedSeq Perl API

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

This package includes a Perl interface to the Genome Analysis Toolkit.

	use NGS::Tools::SpeedSeq;

    my $speedseq = NGS::Tools::SpeedSeq->new(
        bam              => $config,
        ..
    );
    
	my $align = $speedseq->align_seq(); 
    ...

=head1 CLASS/INSTANCE VARIABLES

=cut

has 'config' => ( is => 'ro', );


=head2 $obj->phaseI_root, $obj->phaseII_root

Root directory for storing phase I and phase II output files

=cut

has 'parent_dir' => (
    is  => 'ro',
    isa => 'Str'
);

has 'ALN_root' => (
    is  => 'ro',
    isa => 'Str'
);

has 'VAR_root' => (
    is  => 'ro',
    isa => 'Str'
);

has 'SV_root' => (
    is  => 'ro',
    isa => 'Str'
);

=header2

Various SLURM scheduler jobnames

=cut
has 'jobids' => (
	is => 'ro',
	isa => 'HashRef',
	writer => 'set_jobids',
);

has 'output_files' => (
	is => 'ro',
	isa => 'HashRef',
	writer => 'set_output_files',
);

has 'sample_names' => (
    is     => 'rw',
    isa    => 'ArrayRef',
    writer => 'set_sample_name',
);

=head2 $obj->opts

Attributes related to config

=cut

has 'opts' => (
    is  => 'ro',
    isa => 'HashRef',
);

=head2 $obj->modules

Attibutes for various modules 

=cut

has 'modules' => (
    is  => 'ro',
    isa => 'ArrayRef'
);

has 'dir' => (
    is     => 'ro',
    isa    => 'HashRef',
    writer => 'set_dir_structure',
);

=head1 SUBROUTINES/METHODS
=cut

sub dir_structure {
    my $self = shift;

    my %dirs = (
        ALN => {
            intermediate_bam_files         => 'intermediate_bam',
            merged_bam            => 'merged_bam',
            scripts                        => 'scripts',
            script_alignment     => 'scripts/alignment',
            script_merge_bams => 'scripts/merge_bams',
            log                            => 'log',
            log_alignment        => 'log/alignment',
            log_merge_bams    => 'log/merge_bams',
        },
        SNV => {
            raw_vcf                        => 'raw_vcf',
            scripts                        => 'scripts',
            script_call_SNV => 'scripts/call_SNV',
            log                            => 'log',
            log_call_SNV              => 'log/call_SNV',
        },
        SV => {
            vcf           => 'vcf',
            scripts                     => 'scripts',
            log                         => 'log',
        }
    );

    $self->set_dir_structure( \%dirs );
    $self;
}

=head2 

=cut

sub submit_jobs {
	my $self = shift;
	my $todo = shift;
	my $speedseq_options;
	my $output_files;
	my $jobids;
	
	if ($todo eq 'align'){
		$speedseq_options = $self->get_speedseq_options( step => 'align' );

		# switch to the output directory
		chdir(join( '/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'intermediate_bam_files'}));		

	    foreach my $sample_id ( sort keys %{ $self->config } ) {
			foreach my $lane (keys %{ $self->config->{$sample_id}}){
		
				my $read1 = $self->config->{$sample_id}->{$lane}->{'R1'};
				my $read2 = $self->config->{$sample_id}->{$lane}->{'R2'};
				my $input_string = $read1 . ' '. $read2;
				my $dir = join( '/', $self->opts->{'output_dir'}, $self->dir->{'ALN'}->{'intermediate_bam_files'} );
				my $output_prefix = $sample_id .'_'. $lane;

				# make commandline for slurm script
				my $commandline = $self->make_commandline(
					step        => $todo,
					input       => $input_string,
					output_prefix      => $output_prefix,
					sample_name => $sample_id,
					read_group  => "'\@RG\\tID:$lane\\tSM:$sample_id\\tLB:$sample_id'",
					nct         => $speedseq_options->{'nct'},
					memory      => $speedseq_options->{'mem'}
				);

				my $scriptname = join( '_', 'slurm', 'align_seq', $sample_id, $lane . "\.sh" );
				my $jobname = join( '_', 'alnseq', $sample_id, $lane );
	
				my $dependency = 'none';
				my $slurm_jobid;

				if ( $self->opts->{'align_seq'} eq 'Y' ){
					$slurm_jobid = CCR::Utilities::General::create_submit_slurm_job(
						debug           => $self->opts->{'debug'},
						script          => $scriptname,
						command         => $commandline,
						script_dir      => join( '/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'script_alignment'} ),
						modules         => $self->modules,
						job_name        => $jobname,
						time            => $speedseq_options->{'time'},
						nodes           => 1,
						memory          => $speedseq_options->{'mem'}, 
						ntasks_per_node => $speedseq_options->{'nct'}, 
						dependency      => $dependency,
						output          => join( '/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'log_alignment'}, $jobname )
					);
				} else {
					print "Sequence alignment will not be performed!!\n";
				}
				$output_files->{$sample_id}->{$lane} = join( '/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'intermediate_bam_files'}, $output_prefix . '.bam');
				$jobids->{$todo}->{$sample_id}->{$lane} = $slurm_jobid;
			}
		}
		# switch back to working dir
		chdir($self->parent_dir);

		$self->set_jobids($jobids);
		$self->set_output_files($output_files);
	} elsif($todo eq 'merge_bam') {
		$speedseq_options = $self->get_speedseq_options( step => 'merge_bams' );

		# switch to the output directory
		chdir(join( '/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'merged_bam'}));
		foreach my $sample_id ( sort keys %{ $self->config } ) {
		
	    	foreach my $type ( 'discordants', 'splitters', 'main'){
				my @input_bams;
				foreach my $lane (keys %{ $self->config->{$sample_id}}){
					my $input_bam = join('/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'intermediate_bam_files'}, $sample_id . '_'. $lane .'.'. $type . '.bam');
					$input_bam =~ s/main\.//;
					push @input_bams, $input_bam;
				}
	
				my $commandline = $self->make_commandline(
	                step        => $todo,
	                input       => join(" ", @input_bams),
	                sample_name => $sample_id . '_merged_'. $type,
	                nct         => $speedseq_options->{'nct'},
	                memory      => $speedseq_options->{'mem'}
				);
	
				if ( $self->opts->{'merge_bams'} eq 'Y' ){
					my $scriptname = join( '_', 'slurm', 'merge_bam', $sample_id, $type. "\.sh" );
					my $jobname = join( '_', 'merge_bams', $sample_id, $type );
	
					my $slurm_jobid;
					my $dependency;
					
					if ($self->opts->{'align_seq'} eq 'N') {
						$dependency = 'none';
					}
					elsif (exists $self->jobids->{'align'}->{$sample_id}) {
						my @tmp = values % {$self->jobids->{'align'}->{$sample_id}};
						$dependency = join(':', @tmp);
					}
					
					$slurm_jobid = CCR::Utilities::General::create_submit_slurm_job(
	                    debug           => $self->opts->{'debug'},
	                    script          => $scriptname,
	                    command         => $commandline,
	                    script_dir      => join( '/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'script_merge_bams'} ),
	                    modules         => $self->modules,
	                    job_name        => $jobname,
	                    time            => $speedseq_options->{'time'},
	                    nodes           => 1,
	                    memory          => $speedseq_options->{'mem'},
	                    ntasks_per_node => $speedseq_options->{'nct'},
	                    dependency      => $dependency,
	                    output          => join( '/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'log_merge_bams'}, $jobname )
	                );
					
					my $merged_bam = join('/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'merged_bam'}, $sample_id . '_merged_'. $type .'.bam');
	
					# index the merged bams
					$commandline = join(' ', 'sambamba', 'index', $merged_bam);
					$scriptname = join( '_', 'slurm', 'index_bam', $sample_id, $type . "\.sh" );
					$jobname = join( '_', 'index_bam', $sample_id, $type );
					
					my $slurm_jobid2 = CCR::Utilities::General::create_submit_slurm_job(
	                    debug           => $self->opts->{'debug'},
	                    script          => $scriptname,
	                    command         => $commandline,
	                    script_dir      => join( '/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'script_merge_bams'} ),
	                    modules         => $self->modules,
	                    job_name        => $jobname,
	                    time            => $speedseq_options->{'time'},
	                    nodes           => 1,
	                    memory          => $speedseq_options->{'mem'},
	                    ntasks_per_node => $speedseq_options->{'nct'},
	                    dependency      => $slurm_jobid,
	                    output          => join( '/', $self->parent_dir, 'ALN', $self->dir->{'ALN'}->{'log_merge_bams'}, $jobname )
	                );
	
	
					$output_files->{$todo}->{$sample_id}->{$type} = $merged_bam;
					push @{$jobids->{$todo}->{$sample_id}}, $slurm_jobid2;
				}
			}
		}
		$self->set_jobids($jobids);
		$self->set_output_files($output_files);
		
		# switch back to working dir
		chdir($self->parent_dir);         
	} elsif ($todo eq 'sv'){
		$speedseq_options = $self->get_speedseq_options( step => 'sv' );

		# switch to the output directory
		chdir(join( '/', $self->parent_dir, 'SV', $self->dir->{'SV'}->{'vcf'}));
		foreach my $sample_id ( sort keys %{ $self->config } ) {
			my $input_str = join(' ', '-B', $self->output_files->{'merge_bam'}->{$sample_id}->{'main'}, '-S', $self->output_files->{'merge_bam'}->{$sample_id}->{'splitters'}, '-D', $self->output_files->{'merge_bam'}->{$sample_id}->{'discordants'});
			
			my $commandline = $self->make_commandline(
                step        => $todo,
                input       => $input_str,
                sample_name => $sample_id,
                nct         => $speedseq_options->{'nct'},
                memory      => $speedseq_options->{'mem'}
            );

 			my $scriptname = join( '_', 'slurm', 'sv', $sample_id . "\.sh" );
			my $jobname = join( '_', 'sv', $sample_id );
			my $dependency = 'none';
			if (exists $self->jobids->{'merge_bam'}->{$sample_id}){
				$dependency = join(":", @{$self->jobids->{'merge_bam'}->{$sample_id}});
			}
		
				
			my $slurm_jobid = CCR::Utilities::General::create_submit_slurm_job(
                    debug           => $self->opts->{'debug'},
                    script          => $scriptname,
                    command         => $commandline,
                    script_dir      => join( '/', $self->parent_dir, 'SV', $self->dir->{'SV'}->{'scripts'} ),
                    modules         => $self->modules,
                    job_name        => $jobname,
                    time            => $speedseq_options->{'time'},
                    nodes           => 1,
                    memory          => $speedseq_options->{'mem'},
                    ntasks_per_node => $speedseq_options->{'nct'},
                    dependency      => $dependency,
                    output          => join( '/', $self->parent_dir, 'SV', $self->dir->{'SV'}->{'log'}, $jobname )
            );


			#$output_files->{$todo}->{$sample_id} = join('/', $self->parent_dir, 'SV', $self->dir->{'SV'}->{'vcf'}, $sample_id .'sv.vcf.gz');
			#$jobids->{$todo}->{$sample_id} = $slurm_jobid;
		}
				
		# switch back to working dir
		chdir($self->parent_dir);
	}
	$self;
}

	

sub get_speedseq_options {
    my $self = shift;
    my %args = validate(
        @_,
        {
            step => {
                type     => SCALAR,
                required => 0,
                default  => undef
            }
        }
    );

    my $step   = $args{'step'};

    # initialize with defaul values
    my $time = '72:00:00';
    my $nt   = 1;
    my $nct  = 1;
    my $mem  = 48_000;
    my $rv;
    my $partition = 'general-compute';

        if ( $step eq 'align' ) {
            $nct  = 8;
            $mem  = 48_000;
            $time = '20:00:00';
        }
        elsif ( $step eq 'merge_bams' ) {
            $nct  = 6;
            $mem  = 24_000;
            $time = '10:00:00';
        }
        elsif ( $step eq 'index_bams' ) {
            $nct   = 1;
            $mem  = 4_000;
            $time = '05:00:00';
        }
        elsif ( $step eq 'sv' ) {
            $nct   = 8;
            $mem  = 48_000;
            $time = '72:00:00';
        }
        elsif ( $step eq 'var' ) {
            $time = '20:00:00';
            $nct   = 8;
            $mem  = 48_000;
        }
        elsif ( $step eq 'BaseRecalibrator' ) {
            $nct  = 12;
            $mem  = 48_000;
            $time = '50:00:00';
        }
    $rv->{'nct'}       = $nct;
    $rv->{'mem'}       = $mem;
    $rv->{'time'}      = $time;
    $rv->{'partition'} = $partition;
    return ($rv);
}

sub make_commandline {
    my $self = shift;
    my %args = validate(
        @_,
        {
            memory => {
                type     => SCALAR,
                required => 0,
                default  => 10_000
            },
            step => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            input => {

                #type     => SCALAR || hashref,
                required => 0,
                default  => undef
            },
            output_prefix => {
                required => 1,
                default  => undef
            },
            sample_name => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            nt => {
                type     => SCALAR,
                required => 0,
                default  => 1
            },
            nct => {
                type     => SCALAR,
                required => 0,
                default  => 1
            },
            data_type => {
                type     => SCALAR,
                required => 0,
                default  => 'wes'
            },
			read_group => {
				type => SCALAR,
				required => 0,
				default => undef
			},
        }
    );

    my $params;
    my $step          = $args{'step'};
    my $input         = $args{'input'};
    my $output_prefix        = $args{'output_prefix'};
    my $nt            = $args{'nt'};
    my $nct           = $args{'nct'};
    my $data_type     = $args{'data_type'};
	my $ref = $self->opts->{ref};
	my $read_group = $args{'read_group'};
	my $sample_name = $args{'sample_name'};	
	my $cmd;
	
    if ( $step eq 'align' ) {
        $params = join( ' ', '-o', $output_prefix, '-R', $read_group, '-t', $nct,  $ref, $input);
		$cmd = join( ' ', 'speedseq', $step, $params );
    } elsif ($step eq 'merge_bam'){
		$params = join(' ', $sample_name . '.bam', $input);
		$cmd = join( ' ', 'sambamba merge',  $params );
	} elsif ($step eq 'sv'){
		$params = join(' ' , $input, '-R', $ref, '-g -d -A', '-t', $nct, '-x', $self->opts->{exclude_intervals});
		$cmd = join(' ', 'speedseq', 'sv', $params);
	}
    

    return ($cmd);

}


=head1 AUTHOR

Jianxin Wang, C<< <jw24 at buffalo.edu> >>


=cut

no Moose;

__PACKAGE__->meta->make_immutable;

1;    # End of NGS::Tools::GATK

