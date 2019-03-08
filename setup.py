from setuptools import setup
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "pvp",
    version = "0.1",
    author = "Jerald Yam and Daniel Bogema",
    author_email = "daniel.bogema@dpi.nsw.gov.au",
    description = ("An automated pipleline for generating high-quality variants from parasite sequencing"),
    license = "GPL-3.0",
    keywords = "genomics variant calling parasite",
    url = "https://github.com/jyguojun/T_ori_genome",
    packages=['pvp'],
    scripts=['bin/pvp', 'bin/pvp_index_reference', 'bin/pvp_remove_host_reads', 'bin/pvp_read_alignment', 'bin/pvp_alignment_processing', 'bin/pvp_variant_calling', 'bin/pvp_missingness_filter', 'bin/pvp_rare_allele_filter', 'bin/pvp_create_snp_dict', 'bin/pvp_uniqueness_filter'],
    long_description=read('README.md'),
)
