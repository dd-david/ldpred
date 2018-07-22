#-*- coding: utf-8 -*-
"""
Code for handling plink files.

Uses plinkio.
"""
from plinkio import plinkfile
import scipy as sp


# 유전체(chromosomes) 별 사전을 만듭니다.
# {'chrom_1' : { 'sids' : ['sid_0', 'sid_1', ...],
#                'snp_indices' : [0, 1, ...],
#                'positions' : [0L, 1L, ...],
#                'nts': [['A','G'], ['A','G'], ...]
#              }, ...
# }
def get_chrom_dict(loci, chromosomes):
    chr_dict = {}
    for chrom in chromosomes:
        chr_str = 'chrom_%s' % chrom
        chr_dict[chr_str] = {'sids':[], 'snp_indices':[], 'positions':[], 'nts':[]}

    for i, l in enumerate(loci):
        chrom = l.chromosome
        pos = l.bp_position
        chr_str = 'chrom_%d' % chrom

        chr_dict[chr_str]['sids'].append(l.name)
#         chr_dict[chr_str]['sids'].append('%d_%d'%(chrom,pos))
        chr_dict[chr_str]['snp_indices'].append(i)
        chr_dict[chr_str]['positions'].append(pos)
        chr_dict[chr_str]['nts'].append([l.allele1, l.allele2])

    print 'Genotype dictionary filled'
    return chr_dict

def parse_plink_snps(genotype_file, snp_indices):
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()
    num_individs = len(samples)
    num_snps = len(snp_indices)
    raw_snps = sp.empty((num_snps, num_individs), dtype='int8')
    # If these indices are not in order then we place them in the right place while parsing SNPs.
    snp_order = sp.argsort(snp_indices)
    ordered_snp_indices = list(snp_indices[snp_order])
    ordered_snp_indices.reverse()
    print 'Iterating over file to load SNPs'
    snp_i = 0
    next_i = ordered_snp_indices.pop()
    line_i = 0
    max_i = ordered_snp_indices[0]
    while line_i <= max_i:
        if line_i < next_i:
            plinkf.next()
        elif line_i == next_i:
            line = plinkf.next()
            snp = sp.array(line, dtype='int8')
            bin_counts = line.allele_counts()
            if bin_counts[-1] > 0:
                mode_v = sp.argmax(bin_counts[:2])
                snp[snp == 3] = mode_v
            s_i = snp_order[snp_i]
            raw_snps[s_i] = snp
            if line_i < max_i:
                next_i = ordered_snp_indices.pop()
            snp_i += 1
        line_i += 1
    plinkf.close()
    assert snp_i == len(raw_snps), 'Failed to parse SNPs?'
    num_indivs = len(raw_snps[0])
    freqs = sp.sum(raw_snps, 1, dtype='float32') / (2 * float(num_indivs))
    return raw_snps, freqs

# GWAS Study pipeline
# 링크 : http://www.gwaspi.org/?page_id=123
#
#
# .fam 파일의 형태
# familyId  sampleId  paternalId  maternalId  sex  affection
#       i0        i0          i0          i0    2          1
#       i1        i1          i1          i1    2          2
# sex = [1=male, 2=female, other=unknown]
# affection = [0=unknown, 1=unaffected, 2=affected]
def get_phenotypes(plinkf):

    # 개인(Sample) 정보를 읽어와서, 'phenotype' 이 중복되지 않는 집합을 만듭니다.
    # (ex) LDpred_cc_data_p0.001_test_0 : [0.0, 1.0]
    #      LDpred_data_p0.001_test_0    : [1.3305020332336426, -0.38473600149154663, -1.2226539850234985, ...]
    #
    # sample = [fid, iid, father_iid, mother_iid, sex, affection, phenotype]
    #           ---  ---                          ---  ---------
    #       sex = [0=female, 1=male]
    # affection = [0=control, 1=case, -9=missing] + anything else means that the phenotype is continuous.
    # pehnotype = [0.0=control, 1.0=case]
    samples = plinkf.get_samples()
    num_individs = len(samples)
    Y = [s.phenotype for s in samples]
    fids = [s.fid for s in samples]
    iids = [s.iid for s in samples]
    unique_phens = sp.unique(Y)

    # 'phenotype' 1이면, 형질을 찾을 수 없다고 봐야합니다. (차이가 없음)
    if len(unique_phens) == 1:
        print 'Unable to find phenotype values.'
        has_phenotype = False
    # 'phenotype' 2이면, 차이를 보이는 두 집단을 구분해 낼 수 있습니다.
    # LDpred_cc_data_p0.001_test_0
    elif len(unique_phens) == 2:
        cc_bins = sp.bincount(Y)
        assert len(cc_bins) == 2, 'Problems with loading phenotype'
        print 'Loaded %d controls and %d cases' % (cc_bins[0], cc_bins[1])
        has_phenotype = True
    # 'phenotype' 개수가 다양하게 발견되었다고 볼 수 있습니다.
    # LDpred_data_p0.001_test_0
    else:
        print 'Found quantitative phenotype values'
        has_phenotype = True
    return {'has_phenotype':has_phenotype, 'fids':fids, 'iids':iids, 'phenotypes':Y, 'num_individs':num_individs}
