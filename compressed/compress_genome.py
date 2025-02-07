import sys
from Bio import SeqIO
import csv

def parse_fasta(fasta_file):
    return {record.id: record.seq for record in SeqIO.parse(fasta_file, "fasta")}

def parse_gtf(gtf_file):
    genes = {}
    transcripts = {}
    exons = {}
    with open(gtf_file, 'r') as gtf:
        reader = csv.reader(gtf, delimiter='\t')
        for row in reader:
            if row[0].startswith('#'):
                continue
            attributes = parse_attributes(row[8])
            if row[2] == 'gene':
                if attributes['gene_biotype'] == 'protein_coding':
                    genes[attributes['gene_id']] = {
                        'seqid': row[0],
                        'start': int(row[3]),
                        'end': int(row[4]),
                        'strand': row[6],
                        'gene_name': attributes.get('gene_name', ''),
                        'gene_version': attributes.get('gene_version', ''),
                        'gene_biotype': attributes.get('gene_biotype', '')
                    }
            elif row[2] == 'transcript':
                if attributes['transcript_biotype'] == 'protein_coding' and \
                    (attributes['tag'] == 'Ensembl_canonical' or attributes['tag'] == 'gencode_basic'):
                    transcripts[attributes['transcript_id']] = {
                        'gene_id': attributes['gene_id'],
                        'seqid': row[0],
                        'start': int(row[3]),
                        'end': int(row[4]),
                        'strand': row[6],
                        'gene_name': attributes.get('gene_name', ''),
                        'transcript_biotype': attributes.get('transcript_biotype', ''),
                        'tag' : attributes.get('tag', '')
                    }
            elif row[2] == 'exon':
                if attributes['transcript_biotype'] == 'protein_coding' and \
                    (attributes['tag'] == 'Ensembl_canonical' or attributes['tag'] == 'gencode_basic'):
                    transcript_id = attributes['transcript_id']
                    if transcript_id not in exons:
                        exons[transcript_id] = []
                    exons[transcript_id].append({
                        'start': int(row[3]),
                        'end': int(row[4]),
                        'strand': row[6],
                        'gene_id': attributes['gene_id'],
                        'gene_version': attributes.get('gene_version', ''),
                        'transcript_id': transcript_id,
                        'gene_name': attributes.get('gene_name', ''),
                        'gene_biotype': attributes.get('gene_biotype', ''),
                        'transcript_biotype': attributes.get('transcript_biotype', ''),
                        'exon_id': attributes['exon_id'],
                        'exon_version': attributes.get('exon_version', ''),
                        'tag' : attributes.get('tag', '')
                    })
    return genes, transcripts, exons

def parse_attributes(attribute_string):
    attributes = {}
    for attribute in attribute_string.split(';'):
        if attribute.strip():
            key_value = attribute.strip().split(' ', 1)
            if len(key_value) == 2:
                key, value = key_value
                attributes[key] = value.strip('"')
    return attributes

def get_longest_transcripts(transcripts, exons):
    gene_longest_transcript = {}
    for transcript_id, transcript in transcripts.items():
        gene_id = transcript['gene_id']
        transcript_exons = exons.get(transcript_id, [])
        exonic_length = sum(exon['end'] - exon['start'] for exon in transcript_exons)
        if gene_id not in gene_longest_transcript or exonic_length > gene_longest_transcript[gene_id]['length']:
            gene_longest_transcript[gene_id] = {'transcript_id': transcript_id, 'length': exonic_length}
    return {gene_id: data['transcript_id'] for gene_id, data in gene_longest_transcript.items()}

def create_pseudo_chromosome(gene_longest_transcript, seqs, genes, transcripts, exons):
    # Create a list of transcripts sorted by seq_id and start position
    sorted_transcripts = sorted(
        gene_longest_transcript.items(),
        key=lambda item: (transcripts[item[1]]['seqid'], transcripts[item[1]]['start'])
    )
    pseudo_chrom_seq = {}
    gene_coordinates = {}
    exon_coordinates = {}
    current_position = 1
    seqid_prev = ''
    for gene, transcript_id in sorted_transcripts:
        seqid = transcripts[transcript_id]['seqid']
        if seqid != seqid_prev:
            pseudo_chrom_seq[seqid] = ''
            seqid_prev = seqid
            current_position = 1
        transcript_exons = exons.get(transcript_id, [])
        transcript_exons.sort(key=lambda e: e['start'])
        exon_free_seq = ""
        exon_current_position = current_position
        for exon in transcript_exons:
            exon_id = exon['exon_id']
            exon_new_seq = seqs[seqid][exon['start']-1:exon['end']]
            exon_free_seq += exon_new_seq
            exon_coordinates[exon_id] = (exon_current_position, exon_current_position + len(exon_new_seq) - 1)
            exon_current_position += len(exon_new_seq)
        pseudo_chrom_seq[seqid] += exon_free_seq
        gene_coordinates[gene] = (current_position, current_position + len(exon_free_seq) - 1)
        current_position += len(exon_free_seq)
    return pseudo_chrom_seq, gene_coordinates, exon_coordinates

def write_pseudo_chromosome_fasta(pseudo_chromosomes, output_fasta_file):
    with open(output_fasta_file, "w") as f:
        for chromosome_id, seq in pseudo_chromosomes.items():
            f.write(">{}\n".format(chromosome_id))
            for i in range(0, len(seq), 100):
                f.write(str(seq[i:i+100]) + "\n")

def modify_gff3(genes, transcripts, gene_longest_transcript, pseudo_chromosomes, gene_coordinates, exons, exon_coordinates, output_gff3_file):
    # Create a list of transcripts sorted by seq_id and start position
    sorted_transcripts = sorted(
        gene_longest_transcript.items(),
        key=lambda item: (transcripts[item[1]]['seqid'], transcripts[item[1]]['start'])
    )
    with open(output_gff3_file, "w") as f:
        for gene, transcript_id in sorted_transcripts: #gene_longest_transcript.items():
            gene_start, gene_end = gene_coordinates[gene]
            gene_data = genes[gene]
            f.write("{}\t.\tgene\t{}\t{}\t.\t{}\t.\tID={};gene_name={};gene_version={};gene_biotype={}\n".format(
                transcripts[transcript_id]['seqid'], gene_start, gene_end, gene_data['strand'], gene, gene_data['gene_name'], gene_data['gene_version'], gene_data['gene_biotype']
            ))
            f.write("{}\t.\ttranscript\t{}\t{}\t.\t{}\t.\tID={};Parent={};gene_name={};gene_version={};gene_biotype={};transcript_biotype={};tag={}\n".format(
                transcripts[transcript_id]['seqid'], gene_start, gene_end, gene_data['strand'], transcript_id, gene, gene_data['gene_name'], gene_data['gene_version'], gene_data['gene_biotype'], transcripts[transcript_id]['transcript_biotype'], transcripts[transcript_id]['tag']
            ))
            transcript_exons = exons.get(transcript_id, [])
            transcript_exons.sort(key=lambda e: e['start'])
            for exon in transcript_exons:
                exon_start = exon_coordinates[exon['exon_id']][0]
                exon_end = exon_coordinates[exon['exon_id']][1]
                f.write("{}\t.\texon\t{}\t{}\t.\t{}\t.\tID={};Parent={};gene_id={};gene_name={};gene_version={};gene_biotype={};transcript_biotype={};tag={};exon_version={}\n".format(
                    transcripts[transcript_id]['seqid'], exon_start, exon_end, gene_data['strand'], exon['exon_id'], transcript_id, gene, exon['gene_name'], exon['gene_version'], exon['gene_biotype'], transcripts[transcript_id]['transcript_biotype'], exon['tag'], exon['exon_version']
                ))

def modify_gtf(genes, transcripts, gene_longest_transcript, pseudo_chromosomes, gene_coordinates, exons, exon_coordinates, output_gtf_file):
    # Create a list of transcripts sorted by seq_id and start position
    sorted_transcripts = sorted(
        gene_longest_transcript.items(),
        key=lambda item: (transcripts[item[1]]['seqid'], transcripts[item[1]]['start'])
    )
    with open(output_gtf_file, "w") as f:
        for gene, transcript_id in sorted_transcripts: #gene_longest_transcript.items():
            gene_start, gene_end = gene_coordinates[gene]
            gene_data = genes[gene]
            transcript_data = transcripts[transcript_id]
            gene_name = gene_data.get('gene_name', '')
            gene_version = gene_data.get('gene_version', '')
            gene_biotype = gene_data.get('gene_biotype', '')
            f.write("{}\t.\tgene\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; gene_version \"{}\"; gene_name \"{}\"; gene_biotype \"{}\";\n".format(
                transcript_data['seqid'], gene_start, gene_end, gene_data['strand'], gene, gene_version, gene_name, gene_biotype
            ))
            f.write("{}\t.\ttranscript\t{}\t{}\t.\t{}\t.\ttranscript_id \"{}\"; gene_id \"{}\"; gene_name \"{}\"; gene_version \"{}\"; gene_biotype \"{}\"; transcript_biotype \"{}\"; tag \"{}\"\n".format(
                transcript_data['seqid'], gene_start, gene_end, gene_data['strand'], transcript_id, gene, gene_data['gene_name'], gene_data['gene_version'], gene_data['gene_biotype'], transcripts[transcript_id]['transcript_biotype'], transcripts[transcript_id]['tag']
            ))
            transcript_exons = exons.get(transcript_id, [])
            transcript_exons.sort(key=lambda e: e['start'])
            for exon in transcript_exons:
                exon_start = exon_coordinates[exon['exon_id']][0]
                exon_end = exon_coordinates[exon['exon_id']][1]
                exon_id = exon['exon_id']
                exon_version = exon['exon_version']
                f.write("{}\t.\texon\t{}\t{}\t.\t{}\t.\texon_id \"{}\"; transcript_id \"{}\"; gene_id \"{}\"; gene_name \"{}\"; gene_version \"{}\"; gene_biotype \"{}\"; transcript_biotype \"{}\"; tag \"{}\"; exon_version \"{}\";\n".format(
                    transcript_data['seqid'], exon_start, exon_end, gene_data['strand'], exon_id, transcript_id, gene, gene_name, exon['gene_version'], gene_biotype, transcripts[transcript_id]['transcript_biotype'], exon['tag'], exon_version))

#fasta_file_in, gtf_file_in, fasta_file_out, gtf_file_out, gff_file_out = "mouse_genome.fa", "mouse_genome.gtf", "mouse_genome_compressed.fa", "mouse_genome_compressed.gtf", "mouse_genome_compressed.gff3"
fasta_file_in, gtf_file_in, fasta_file_out, gtf_file_out, gff_file_out = sys.argv[1:]

print("loading seqs")
seqs = parse_fasta(fasta_file_in)
print("parsing gtf")
genes, transcripts, exons = parse_gtf(gtf_file_in)
print("finding longest transcript")
gene_longest_transcript = get_longest_transcripts(transcripts, exons)
print("creating pseudo-chromosomes")
pseudo_chromosomes, gene_coordinates, exon_coordinates = create_pseudo_chromosome(gene_longest_transcript, seqs, genes, transcripts, exons)
print("writing output files")
write_pseudo_chromosome_fasta(pseudo_chromosomes, fasta_file_out) #, chromosome_id)
modify_gff3(genes, transcripts, gene_longest_transcript, pseudo_chromosomes, gene_coordinates, exons, exon_coordinates, gff_file_out)
modify_gtf(genes, transcripts, gene_longest_transcript, pseudo_chromosomes, gene_coordinates, exons, exon_coordinates, gtf_file_out)

