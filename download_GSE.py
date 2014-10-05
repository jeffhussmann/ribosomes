import argparse
import ftplib
import urlparse
import subprocess
import os
import xml.etree.ElementTree as etree
import string

def get_xml(paper_dir, accession):
    '''Download an xml file containing information about a GEO accession number.
    '''
    xml_tail = '{0}_family.xml'.format(accession)
    xml_fn = '{0}/{1}'.format(paper_dir, xml_tail)
    tgz_fn = xml_fn + '.tgz'

    xml_url_template = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{accession}/miniml/{xml_tail}.tgz'
    xml_url = xml_url_template.format(prefix=accession[:5],
                                      accession=accession,
                                      xml_tail=xml_tail,
                                     )

    if not os.path.isdir(paper_dir):
        os.makedirs(paper_dir)

    wget_command = ['wget', '--quiet', '-P', paper_dir, xml_url]
    subprocess.check_call(wget_command)

    tar_command = ['tar', 'xzf', tgz_fn, '-C', paper_dir]
    subprocess.check_call(tar_command)
    os.remove(tgz_fn)
    os.remove(xml_fn)

    return xml_fn

def extract_samples_from_xml(xml_fn, condition=lambda x: True):
    '''Parse an xml file describing a GEO accession number to extract samples.
    
    Sanitizes the names of samples by replacing spaces, slashes, and brackets
    with underscores, replacing unicode deltas with the string 'delta_', and 
    removing parentheses.

    condition -- a filtering function to be applied to sample names
    
    Returns a list of (sample name, list of URLs) tuples.
    '''
    samples = []
    tree = etree.parse(xml_fn)
    root = tree.getroot()
    sanitize_table = string.maketrans(' /,[]', '_____')
    for child in root:
        if child.tag.endswith('Sample'):
            skip_sample = False
            for grand in child:
                if grand.tag.endswith('Title'):
                    sample_name = grand.text
                    sample_URLs = []
                    if not type(sample_name) == str:
                        # Some sample names (knockout strains) have a unicode
                        # delta in them. Replace this with a spelled out
                        # 'delta' and coerce to a non-unicode string.
                        delta = u'\u0394'
                        sample_name = str(sample_name.replace(delta, 'delta_'))
                    sample_name = sample_name.translate(sanitize_table, '()')
                    if not condition(sample_name):
                        skip_sample = True
                        break
                elif grand.tag.endswith('Supplementary-Data') and grand.attrib['type'] == 'SRA Experiment':
                    sample_URL = grand.text.strip()
                    sample_URLs.append(sample_URL)
            if not skip_sample:
                samples.append((sample_name, sample_URLs))
    return samples

def get_run_info(run, paper_dir):
    url = 'http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?retmode=xml&run={}'.format(run)
    xml_fn = '{0}/{1}.xml'.format(paper_dir, run)

    wget_command = ['wget', '--quiet', url, '-O', xml_fn]
    subprocess.check_call(wget_command)

    tree = etree.parse(xml_fn)
    root = tree.getroot()

    info = {}
    info['size'] = int(root.find('RUN').attrib['size'])

    layout = root.find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_LAYOUT')
    if layout.find('SINGLE') is not None:
        info['layout'] = 'single'
    elif layout.find('PAIRED') is not None:
        info['layout'] = 'paired'
        info['nominal_length'] = layout.find('PAIRED').get('NOMINAL_LENGTH')
    else:
        info['layout'] = 'unknown'
    
    os.remove(xml_fn)
    return info

def download_samples(paper_dir, samples):
    '''Downloads samples using ascp.

    Requires the environment variable ASCP_KEY to be set to a path to the key.
    '''
    sra_fns = []
    for sample_name, sample_URLs in samples:
        sample_dir = '{0}/{1}/data'.format(paper_dir, sample_name)
        if not os.path.isdir(sample_dir):
            os.makedirs(sample_dir)
        
        for sample_URL in sample_URLs:
            parsed_sample_url = urlparse.urlparse(sample_URL)
            f = ftplib.FTP(parsed_sample_url.netloc)
            f.login()
            f.cwd(parsed_sample_url.path)
            ls = []
            f.dir(ls.append)
            f.quit()
            
            for line in ls:
                run = line.split()[-1]
                info = get_run_info(run, paper_dir)
                if info['size'] > 15e9:
                    continue
                else:
                    print run, info['size']
                run_url = '{0}/{1}/{1}.sra'.format(sample_URL, run)
                parsed_run_url = urlparse.urlparse(run_url)
                wget_command = ['wget',
                                '-P', sample_dir,
                                run_url,
                               ]
                ascp_command = ['ascp',
                                '-i', os.environ['ASCP_KEY'],
                                '-QT',
                                '-l', '300m',
                                'anonftp@{0}:{1}'.format(parsed_run_url.netloc, parsed_run_url.path),
                                sample_dir,
                               ]
                subprocess.check_call(ascp_command)
                sra_fn = '{0}/{1}.sra'.format(sample_dir, run)
                sra_fns.append((sra_fn, info['layout']))

    return sra_fns

def dump_fastqs(sra_fns):
    '''Dumps fastq files from sra files, then deletes the sra files. ''' 
    for sra_fn, layout in sra_fns:
        print "fastq-dump'ing {0}".format(sra_fn) 
        head, tail = os.path.split(sra_fn)
        root, ext = os.path.splitext(sra_fn)
        if layout == 'paired':
            # Split into two files and include read number (out of pair) in the
            # seq name line.
            fastq_dump_command = ['fastq-dump',
                                  '--split-3',
                                  '--defline-seq', '@$ac.$si.$ri',
                                  '--defline-qual', '+',
                                  '-O', head,
                                  sra_fn,
                                 ]
        elif layout == 'single':
            fastq_dump_command = ['fastq-dump',
                                  '--defline-seq', '@$ac.$si',
                                  '--defline-qual', '+',
                                  '-O', head,
                                  sra_fn,
                                 ]
        else:
            raise ValueError('layout not known')
        
        subprocess.check_call(fastq_dump_command)
        os.remove(sra_fn)

experiments = {
    # Yeast
    # Ribosome profiling
    'guydosh_cell': ('GSE52968', lambda name: 'wild' not in name and 'disome' in name),
    'brar_science': ('GSE34082', lambda name: 'mRNA' in name and 'exponential' in name),
    'artieri': ('GSE50049', lambda name: 'Mixed_parental' in name),
    'mcmanus_gr': ('GSE52119', lambda name: 'S._cerevisiae' in name),
    'ingolia_science': ('GSE13750', lambda name: True),
    'zinshteyn_plos_genetics': ('GSE45366', lambda name: 'WT' in name),
    'lareau_elife': ('GSE58321', lambda name: True),
    'baudin-baillieu_cell_reports': ('GSE41590', lambda name: True),
    'gerashchenko_nar': ('GSE59573', lambda name: 'unstressed' in name),
    # UTR boundary identification
    'arribere_gr': ('GSE39074', lambda name: 'S288C' in name and '_TLSeq' in name),
    'pelechano_nature': ('GSE39128', lambda name: name.startswith('nypd')),
    'three_p_seq': ('GSE53310', lambda name: True),
    'park_nar': ('GSE49026', lambda name: True),
    'wilkening_nar': ('GSE40110', lambda name: 'rnaseq' in name and 'ypd' in name),
    # RNA-seq
    'nagalakshmi_science': ('GSE11209', lambda name: True),
    'yassour_gb': ('GSE21739', lambda name: name == 'dUTP'),
    
    # Mouse
    'guo_nature': ('GSE22004', lambda name: 'wild-type_runs' in name),
    'ingolia_cell': ('GSE30839', lambda name: '_chx' in name),
    'thoreen_nature': ('GSE36892', lambda name: 'Footprint_for_vehicle-treated_wild-type_MEFs_Replicate_1' in name),

    # E. coli
    # barcoding for counting
    'shiroguchi': ('GSE34449', lambda name: True),
    # ribosome profiling
    'li_nature': ('GSE35641', lambda name: 'E._coli_MG1655' in name),
    'li_cell': ('GSE53767', lambda name: True),
}

non_GSE_experiments = {
    'geisler_mc': [('WT_plus', ['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX106/SRX106067']),
                  ],
    'CHM1htert': [('resequencing', ['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP017/SRP017546']),
                 ],
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('paper_name')
    parser.add_argument('papers_dir', help='base directory')
    parser.add_argument('--paired', help='data is paired-end reads', action='store_true')
    parser.add_argument('--list', help='only list samples, don\'t download', action='store_true')
    args = parser.parse_args()
    paper_dir = '{0}/{1}'.format(args.papers_dir, args.paper_name)

    if args.paper_name in experiments:
        accession, condition = experiments[args.paper_name]
        xml_fn = get_xml(paper_dir, accession)
        samples = extract_samples_from_xml(xml_fn, condition=condition)
    elif args.paper_name in non_GSE_experiments:
        samples = non_GSE_experiments[args.paper_name]

    for sample in samples:
        print sample[0]

    if not args.list:
        sra_fns = download_samples(paper_dir, samples)
        dump_fastqs(sra_fns)
