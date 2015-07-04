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
                        deltas = [u'\u0394', u'\u2206']
                        for delta in deltas:
                            sample_name = sample_name.replace(delta, 'delta_')
                        sample_name = str(sample_name)
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

def get_run_info(run, paper_dir, verbose=False):
    url = 'http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?retmode=xml&run={}'.format(run)
    xml_fn = '{0}/{1}.xml'.format(paper_dir, run)

    wget_command = ['wget', '--quiet', url, '-O', xml_fn]
    subprocess.check_call(wget_command)

    tree = etree.parse(xml_fn)
    root = tree.getroot()

    info = {}
    info['run'] = run
    info['size'] = int(root.find('RUN').attrib['size'])
    info['total_spots'] = int(root.find('RUN').attrib['total_spots'])

    layout = root.find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_LAYOUT')
    if layout.find('SINGLE') is not None:
        info['layout'] = 'single'
    elif layout.find('PAIRED') is not None:
        info['layout'] = 'paired'
        info['nominal_length'] = int(layout.find('PAIRED').get('NOMINAL_LENGTH', 0))
    else:
        info['layout'] = 'unknown'
    
    if verbose:
        print 'Run info for {0}:'.format(run)
        print '\tlayout = {0}'.format(info['layout'])
        if 'nominal_length' in info:
            print '\tnominal length = {0:,}'.format(info['nominal_length'])
        print '\tsize = {0:,}'.format(info['size'])
        print '\ttotal spots = {0:,}'.format(info['total_spots'])

    os.remove(xml_fn)
    return info

def download_samples(paper_dir, samples, condition=lambda x: True):
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
                sra_id, _ = os.path.splitext(run)
                info = get_run_info(run, paper_dir, verbose=True)
                if not condition(info):
                    continue
                run_url = '{0}/{1}/{1}.sra'.format(sample_URL, sra_id)
                parsed_run_url = urlparse.urlparse(run_url)
                ascp_command = ['ascp',
                                '-i', os.environ['ASCP_KEY'],
                                '-QT',
                                '-l', '300m',
                                'anonftp@{0}:{1}'.format(parsed_run_url.netloc, parsed_run_url.path),
                                sample_dir,
                               ]
                subprocess.check_call(ascp_command)
                sra_fn = '{0}/{1}.sra'.format(sample_dir, sra_id)
                sra_fns.append((sra_fn, info['layout']))

    return sra_fns

def dump_fastqs(sra_fns, full_names=False, gzip=False):
    '''Dumps fastq files from sra files, then deletes the sra files.
    full_names controls which --defline-seq value is used.
    ''' 
    for sra_fn, layout in sra_fns:
        print "fastq-dump'ing {0}".format(sra_fn) 
        head, tail = os.path.split(sra_fn)
        root, ext = os.path.splitext(sra_fn)
        fastq_dump_command = ['fastq-dump']
        if gzip:
            fastq_dump_command.append('--gzip')
        if layout == 'paired':
            # Split into two files and include read number (out of pair) in the
            # seq name line.
            if full_names:
                name_format = '@$sn.$ri'
            else:
                name_format = '@$ac.$si.$ri'
            fastq_dump_command.extend(['--split-3',
                                       '--defline-seq', name_format,
                                      ])
        elif layout == 'single':
            if full_names:
                name_format = '@$sn'
            else:
                name_format = '@$ac.$si'
            fastq_dump_command.extend(['--defline-seq', name_format])
        else:
            raise ValueError('layout not known')

        fastq_dump_command.extend(['--defline-qual', '+',
                                   '--outdir', head,
                                   sra_fn,
                                  ])
        
        subprocess.check_call(fastq_dump_command)
        os.remove(sra_fn)

experiments = {
    # Yeast
    # Ribosome profiling
    'guydosh_cell': ('GSE52968', lambda name: True),
    'brar_science': ('GSE34082', lambda name: 'exponential' in name and 'traditional' in name),
    'artieri': ('GSE50049', lambda name: 'Mixed_parental' in name),
    'mcmanus_gr': ('GSE52119', lambda name: 'S._cerevisiae' in name),
    'ingolia_science': ('GSE13750', lambda name: True),
    'zinshteyn_plos_genetics': ('GSE45366', lambda name: True),
    'lareau_elife': ('GSE58321', lambda name: True),
    'baudin-baillieu_cell_reports': ('GSE41590', lambda name: True),
    'gerashchenko_nar': ('GSE59573', lambda name: 'oxidative' not in name and 'unstressed' not in name and 'heatshock' not in name),
    'pop_msb': ('GSE63789', lambda name: True),
    'gardin_elife': ('GSE51164', lambda name: True),
    'nedialkova_cell': ('GSE67387', lambda name: 'tEKQ' in name or 'empty_plasmid' in name),
    'jan_science': ('GSE61012', lambda name: name.startswith('BirAmVenus') or name.startswith('BirAHeh')),
    'williams_science': ('GSE61011', lambda name: True),
    'sen_gr': ('GSE66411', lambda name: 'ribo_wild' in name),
    # UTR boundary identification
    'arribere_gr': ('GSE39074', lambda name: 'S288C' in name and '_TLSeq' in name),
    'pelechano_nature': ('GSE39128', lambda name: name.startswith('nypd')),
    'three_p_seq': ('GSE53310', lambda name: True),
    'park_nar': ('GSE49026', lambda name: True),
    'wilkening_nar': ('GSE40110', lambda name: 'rnaseq' in name and 'ypd' in name),
    # RNA-seq
    'nagalakshmi_science': ('GSE11209', lambda name: True),
    'yassour_gb': ('GSE21739', lambda name: name == 'dUTP'),
    # Other
    'pelechano_cell': ('GSE63120', lambda name: 'by_chx_1' in name),
    
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
    'artieri_gr_2': [('non_multiplexed', ['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP033/SRP033885']),
                    ],
    'gardin_elife': [('Sc-His_mRNA-seq', ['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX645/SRX645892']),
                     ('Sc-His_Footprint-seq', ['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX645/SRX645893']),
                     ('Sc-Lys_mRNA-seq', ['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX645/SRX645894']),
                     ('Sc-Lys_Footprint-seq', ['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX645/SRX645895']),
                    ],
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('paper_name')
    parser.add_argument('papers_dir', help='base directory')
    parser.add_argument('--list', help='only list samples, don\'t download', action='store_true')
    parser.add_argument('--full_names', help='use original read names', action='store_true')
    parser.add_argument('--gzip', help='gzip fastqs', action='store_true')
    args = parser.parse_args()
    paper_dir = '{0}/{1}'.format(args.papers_dir.rstrip('/'), args.paper_name)

    samples = []
    if args.paper_name in experiments:
        accession, condition = experiments[args.paper_name]
        xml_fn = get_xml(paper_dir, accession)
        samples.extend(extract_samples_from_xml(xml_fn, condition=condition))
        os.remove(xml_fn)
    if args.paper_name in non_GSE_experiments:
        samples.extend(non_GSE_experiments[args.paper_name])
    if not samples:
        raise ValueError('No samples found for {0}'.format(args.paper_name))

    for sample in sorted(samples):
        print sample[0]

    if not args.list:
        sra_fns = download_samples(paper_dir, samples)
        dump_fastqs(sra_fns, gzip=args.gzip, full_names=args.full_names)
