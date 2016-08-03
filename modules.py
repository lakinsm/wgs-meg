#!/usr/bin/env python3
"""MEG Pipeline for Whole Genome Sequencing"""
from jip import *


@tool(inputs=['input', 'basename'])
def trimmomatic(tool):
    """\
    Quality filtering on input reads with Trimmomatic

    usage:
        trimmomatic [-P <pipeline>] [-T <threads>] -i <inputs>... -b <basename>

    Options:
        -P <pipeline>			SE or PE for single end and paired end, respectively [default: PE]
        -T <threads>			Number of threads to use [default: 1]
        -i, --input <inputs>		List of input .fastq.gz or .fastq files
        -b, --basename <basename>	Basename for output files
    """
    def validate(self):
        if tool.pipeline == "PE":
            self.add_output('output_forward_paired', self.basename + '_1P')
            self.add_output('output_reverse_paired', self.basename + '_2P')
            self.add_output('output_forward_unpaired', self.basename + '_1U')
            self.add_output('output_reverse_unpaired', self.basename + '_2U')
        elif tool.pipeline == "SE":
            self.add_output('output', self.basename + '_trimmed')

    
    if tool.pipeline == "PE":
        return '''\
        java -jar /usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar \
        ${pipeline} \
        -threads ${threads} \
        %s %s \
        -baseout ${basename} \
        ILLUMINACLIP:/usr/local/Trimmomatic-0.36/adapters/TruSeq3-${pipeline}.fa:2:30:10:3:TRUE \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        '''.format(self.inputs[0], self.inputs[1])
    elif tool.pipeline == "SE":
        return '''\
        java -jar /usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar \
        ${pipeline} \
        -threads ${threads} \
        ${inputs} \
        ${basename}_trimmed \
        ILLUMINACLIP:/usr/local/Trimmomatic-0.36/adapters/TruSeq3-${pipeline}.fa:2:30:10:3:TRUE \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        '''


@pipeline('assembly')
class WGSMegPipeline(object):
    """\
    Run the WGS MEG Pipeline for assembly and alignment characterization.

    usage:
        assembly [-T <threads>] [-P <pipeline>] -i <inputs>... -s <sample_name> -o <output>

    Options:
        -P <pipeline>			SE or PE for single end and paired end, respectively [default: PE]
        -T <threads>            	Number of threads to use [default: 1]
        -i, --input <inputs>		List of input .fastq.gz or .fastq files
        -s, --sample <sample_name>	Basename for sample
        -o, --output <output>		Basename for output directory
    """

    def pipeline(self):
        p = Pipeline()
        ref = p.run('trimmomatic',
                    pipeline=self.P,
                    threads=self.T,
                    inputs=self.input,
                    output=self.sample)
        p.context(locals())
        return p

