import pysam

def update_gedmatch_with_vcf(vcf_file, template_file, output_file, log_file):
    """
    Update Gedmatch template with actual genotypes from VCF
    
    Template contains all array positions with REF/REF calls
    VCF contains variants - we update the template where variants exist
    """
    
    # Statistics counters
    stats = {
        'total_template_positions': 0,
        'found_in_vcf': 0,
        'ref_calls_kept': 0,
        'het_calls': 0,
        'hom_alt_calls': 0,
        'hom_ref_calls': 0,
        'missing_genotypes': 0,
        'non_snp_skipped': 0,
        'updated': 0
    }
    
    with open(log_file, 'w') as log:
        log.write("Starting Gedmatch template update\n")
        
        # Load VCF into memory indexed by rsID
        log.write("Loading VCF...\n")
        vcf_data = {}
        with pysam.VariantFile(vcf_file) as vcf:
            for record in vcf:
                if record.id and record.id.startswith('rs'):
                    vcf_data[record.id] = record
        
        stats['found_in_vcf'] = len(vcf_data)
        log.write(f"Loaded {len(vcf_data)} rsIDs from VCF\n")
        
        # Process template and update with VCF data
        log.write("Processing template...\n")
        with open(template_file, 'r') as template, \
             open(output_file, 'w') as out:
            
            # Copy header
            header = template.readline()
            out.write(header)
            stats['total_template_positions'] = 0
            
            for line in template:
                stats['total_template_positions'] += 1
                parts = line.strip().split('\t')
                
                if len(parts) < 5:
                    continue
                
                rsid, chrom, pos, allele1, allele2 = parts[0], parts[1], parts[2], parts[3], parts[4]
                
                # Check if this rsID has a variant in the VCF
                if rsid in vcf_data:
                    record = vcf_data[rsid]
                    
                    # Get genotype from first sample
                    sample = list(record.samples.keys())[0]
                    gt = record.samples[sample]['GT']
                    
                    # Skip missing genotypes - keep template REF call
                    if None in gt or gt == (None, None):
                        stats['missing_genotypes'] += 1
                        out.write(line)
                        continue
                    
                    # Get alleles
                    alleles = [record.alleles[g] for g in gt]
                    new_allele1 = alleles[0] if len(alleles) > 0 else allele1
                    new_allele2 = alleles[1] if len(alleles) > 1 else new_allele1
                    
                    # Skip if not a SNP (single base) - keep template
                    if len(new_allele1) != 1 or len(new_allele2) != 1:
                        stats['non_snp_skipped'] += 1
                        out.write(line)
                        continue
                    
                    # Count genotype types
                    if gt == (0, 0):
                        stats['hom_ref_calls'] += 1
                        # Keep template REF call
                        out.write(line)
                    elif gt[0] != gt[1]:  # Heterozygous
                        stats['het_calls'] += 1
                        stats['updated'] += 1
                        out.write(f"{rsid}\t{chrom}\t{pos}\t{new_allele1}\t{new_allele2}\n")
                    else:  # Homozygous ALT
                        stats['hom_alt_calls'] += 1
                        stats['updated'] += 1
                        out.write(f"{rsid}\t{chrom}\t{pos}\t{new_allele1}\t{new_allele2}\n")
                else:
                    # Not in VCF - keep template REF call
                    stats['ref_calls_kept'] += 1
                    out.write(line)
        
        # Write statistics to log
        log.write("\n=== Conversion Statistics ===\n")
        log.write(f"Total template positions: {stats['total_template_positions']}\n")
        log.write(f"rsIDs found in VCF: {stats['found_in_vcf']}\n")
        log.write(f"Positions updated: {stats['updated']}\n")
        log.write(f"\nGenotype breakdown:\n")
        log.write(f"  - REF calls kept from template: {stats['ref_calls_kept']}\n")
        log.write(f"  - HOM REF in VCF (0/0): {stats['hom_ref_calls']}\n")
        log.write(f"  - HET calls updated (0/1, 1/0): {stats['het_calls']}\n")
        log.write(f"  - HOM ALT calls updated (1/1): {stats['hom_alt_calls']}\n")
        log.write(f"\nSkipped:\n")
        log.write(f"  - Missing genotypes: {stats['missing_genotypes']}\n")
        log.write(f"  - Non-SNPs (indels): {stats['non_snp_skipped']}\n")
        
        update_rate = (stats['updated'] / stats['total_template_positions'] * 100) if stats['total_template_positions'] else 0
        log.write(f"\nUpdate rate: {update_rate:.2f}% of template positions modified\n")

if __name__ == "__main__":
    update_gedmatch_with_vcf(
        snakemake.input.vcf,
        snakemake.input.template,
        snakemake.output[0],
        snakemake.log[0]
    )