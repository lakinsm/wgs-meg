#!/usr/bin/env python3
"""MEG Pipeline for Whole Genome Sequencing"""
from jip import *


@tool('trimmomatic')
class Trimmomatic(object):
    """\
    Quality filtering on input reads with Trimmomatic

    usage:
        trimmomatic [-P <pipeline>] [-T <threads>] [-i <inputs>...] [-o <output>]

    Options:
        -P <pipeline>       SE or PE for single end and paired end, respectively [default: PE]
        -T <threads>        Number of threads to use [default: 1]
        -i, --input <inputs>     List of input .fastq.gz or .fastq files
        -o, --output <output>    Basename for output files
    """
    def validate(self):
        self.add_output('output_forward_paired', self.output + '_1P')
        self.add_output('output_forward_paired', self.output + '_2P')
        self.add_output('output_forward_paired', self.output + '_1U')
        self.add_output('output_forward_paired', self.output + '_2U')

    def get_command(self):
        return '''\
        java -jar /usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar \
        ${pipeline} \
        -threads ${threads} \
         -baseout ${output} \
         ILLUMINACLIP:/usr/local/Trimmomatic-0.36/adapters/TruSeq3-${pipeline}.fa:2:30:10:3:TRUE \
         LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        '''


@pipeline('wgs-meg')
class WGSMegPipeline(object):
    """\
    Run the WGS MEG Pipeline for assembly and alignment characterization.

    usage:
        wgs-meg [-T <threads>] [-P <pipeline>] -input <inputs>... -output <output>

    Options:
        -P <pipeline>           SE or PE for single end and paired end, respectively [default: PE]
        -T <threads>            Number of threads to use [default: 1]
        -input <inputs>         List of input .fastq.gz or .fastq files
        -sample <sample_name>   Basename for sample
        -output <output>        Basename for output directory
    """

    def pipeline(self):
        p = Pipeline()
        ref = p.run('trimmomatic',
                    pipeline=self.pipeline,
                    threads=self.threads,
                    inputs=self.inputs,
                    output=self.sample_name)
        p.context(locals())
        return p

