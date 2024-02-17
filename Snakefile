
genes = ["pro", "RT-p66", "RT-p51", "INT",
         "p17", "p24", "p2", "p7", "p1", "p6",
         "vif", "vpr", "tat",  "rev", "vpu",
         "gp120", "gp41", "nef"]


rule assemble:
    input:
        "config/reference.fasta",
        "config/genome_annotation.gff3",
        "config/pathogen.json",
        "config/CHANGELOG.md",
        "config/README.md",
        "results/tree.json",
        "results/example_sequences.fasta"
    output:
        directory("dataset")
    shell:
        """
        mkdir -p {output}
        cp {input} {output}
        mv dataset/example_sequences.fasta dataset/sequences.fasta
        """

subtypes_to_exclude = ['U', 'UO', 'O','N', 'P', 'CPZ', 'GOR']

rule name_by_accession:
    message: """renaming sequences by accession"""
    input:
        "data/HIV1_SFL_2021_genome_DNA.fasta"
    output:
        "data/sequences_renamed.fasta"
    run:
        from Bio import SeqIO
        accessions = set()
        with open(output[0], "w") as f:
            for fname in input:
                for record in SeqIO.parse(fname, "fasta"):
                    subtype = record.id.split(".")[0]
                    record.id = record.id.split(".")[-1]
                    if record.id in accessions:
                        continue
                    if subtype in subtypes_to_exclude:
                        continue
                    accessions.add(record.id)
                    record.description = ""
                    record.seq = record.seq.replace('-', '')
                    SeqIO.write(record, f, "fasta")

rule align:
    message: """aligning sequences to the reference using nextclade v3"""
    input:
        sequences = "data/sequences_renamed.fasta",
        ref = "config/reference.fasta",
        annotation = "config/genome_annotation.gff3",
        pjson = "config/pathogen.json"
    output:
        alignment = "results/alignment.fasta",
        translations = expand("results/translation_{gene}.fasta", gene=genes)
    params:
        translations = lambda w:"results/translation_{cds}.fasta"
    threads: 4
    shell:
        """
        nextclade3 run -j {threads} --input-ref {input.ref} \
                   --input-pathogen-json {input.pjson} --input-annotation {input.annotation} \
                   --output-fasta {output.alignment} --output-translations {params.translations} \
                   --include-reference --output-csv results/nextclade.csv  \
                   {input.sequences}
        """

rule mask_for_tree:
    input:
        sequences = "results/alignment.fasta",
        mask = "config/mask_for_tree.bed"
    output:
        masked = "results/masked.fasta"
    shell:
        """
        augur mask --mask {input.mask} --sequences {input.sequences} --output {output.masked}
        """

rule make_metadata:
    message: """pull metadata from description of sequences"""
    input:
        "data/HIV1_SFL_2021_genome_DNA.fasta"
    output:
        "results/metadata.tsv"
    run:
        from Bio import SeqIO
        import pandas as pd
        metadata = [{'strain':'NC_001802','subtype':'B', 'country':'FR', 'date':1983, 'accession':'NC_001802', 'name':"HXB2_reference"}]
        fields = {'strain': -1, 'subtype': 0, 'country': 1, 'date': 2, 'accession': -1,  'name': [3, -1]}
        accessions = {metadata[0]['strain']}
        for fname in input:
            for record in SeqIO.parse(fname, "fasta"):
                entries = record.id.split(".")
                if len(entries)<5:
                    print(entries)
                    continue
                datum = {k: entries[v] if type(v)==int else '.'.join(entries[v[0]:v[1]]) for k, v in fields.items()}
                if datum['subtype'] in subtypes_to_exclude:
                    continue
                try:
                    if datum['date']< '30':
                        datum['date'] = 2000 + int(datum['date'])
                    else:
                        datum['date'] = 1900 + int(datum['date'])
                except:
                    datum['date'] = None
                if datum['strain'] in accessions:
                    continue
                accessions.add(datum['strain'])
                metadata.append(datum)
        pd.DataFrame(metadata).to_csv(output[0], sep='\t', index=False)

max_count = 50
rule split_by_subtype:
    input:
        alignment = "results/masked.fasta",
        metadata = "results/metadata.tsv"
    output:
        alignments = directory("results/subtype_alignments"),
        touch_file = 'results/subtype_alignments_touch_file'
    run:
        import pandas as pd
        import os
        from Bio import SeqIO
        os.makedirs(output.alignments)
        from collections import defaultdict
        alignments_by_subtype = defaultdict(list)
        metadata = pd.read_csv(input.metadata, sep='\t', index_col='strain')
        for record in SeqIO.parse(input.alignment, "fasta"):
            try:
                subtype = metadata.loc[record.id,'subtype']
                alignments_by_subtype[subtype].append(record)
            except KeyError:
                print(record.id, 'not found')

        from random import shuffle
        with open(output.examples, "w") as f_examples:
            for subtype, seqs in alignments_by_subtype.items():
                if len(seqs)<3: continue
                shuffle(seqs)
                with open(output.alignments + f"/{subtype}.fasta", "w") as f:
                    if len(seqs)<max_count:
                        SeqIO.write(seqs, f, "fasta")
                        SeqIO.write(seqs[0], f_examples, "fasta")
                    else:
                        SeqIO.write(seqs[:max_count], f, "fasta")
                        SeqIO.write(seqs[max_count:max_count+3], f_examples, "fasta")

        os.system('touch '+output.touch_file)


rule make_examples:
    input:
        sequences = "data/sequences_renamed.fasta",
        metadata = "results/metadata.tsv",
        tree = "results/tree.nwk",
    output:
        examples = 'results/example_sequences.fasta'
    run:
        import pandas as pd
        import os
        from Bio import SeqIO, Phylo
        from collections import defaultdict
        sequences_by_subtype = defaultdict(list)
        tips = set([n.name for n in Phylo.read(input.tree, "newick").get_terminals()])
        metadata = pd.read_csv(input.metadata, sep='\t', index_col='strain')
        for record in SeqIO.parse(input.sequences, "fasta"):
            if record.id in tips:
                continue
            try:
                mdata = metadata.loc[record.id].to_dict()
                subtype = mdata['subtype']
                record.id = f"{subtype}.{mdata['country']}.{mdata['date']}.{mdata['name']}.{record.id}"
                sequences_by_subtype[subtype].append(record)
            except KeyError:
                print(record.id, 'not found')

        from random import shuffle
        with open(output.examples, "w") as f_examples:
            for subtype, seqs in sequences_by_subtype.items():
                if len(seqs)<2: continue
                shuffle(seqs)
                if len(seqs)<max_count:
                    SeqIO.write(seqs[0], f_examples, "fasta")
                else:
                    SeqIO.write(seqs[:5], f_examples, "fasta")


rule trees_by_subtype:
    input:
        touch_file = "results/subtype_alignments_touch_file"
    output:
        trees = directory("results/subtype_trees"),
        touch_file = "results/subtype_trees_touch_file"
    params:
        alignments = "results/subtype_alignments"
    shell:
        """
        mkdir -p {output.trees}
        for f in {params.alignments}/*.fasta; do
            augur tree --alignment $f --output {output.trees}/$(basename $f .fasta).nwk
        done
        touch {output.touch_file}
        """

rule group_trees:
    message: """collect all trees and attach them as children to one large polytomy"""
    input:
        touch_file = "results/subtype_trees_touch_file"
    output:
        "results/tree_raw.nwk"
    params:
        trees = "results/subtype_trees"
    run:
        from Bio import Phylo
        import glob


        T = Phylo.BaseTree.Tree()
        for t in glob.glob(params.trees + "/*.nwk"):
            sub_tree = Phylo.read(t, "newick")
            sub_tree.root_at_midpoint()
            sub_tree.root.branch_length = 0.1
            T.root.clades.append(sub_tree.root)
        Phylo.write(T, output[0], "newick")


rule refine:
    input:
        tree = "results/tree_raw.nwk",
        metadata = "results/metadata.tsv"
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    shell:
        """
        augur refine --tree {input.tree} --metadata {input.metadata} \
                     --keep-root \
                     --output-tree {output.tree} --output-node-data {output.node_data}
        """

rule clades:
    input:
        tree = "results/tree.nwk",
        metadata = "results/metadata.tsv"
    output:
        "results/clades.json"
    shell:
        """
        python3 scripts/assign_clades.py --tree {input.tree} --metadata {input.metadata} --output {output}
        """


rule ancestral:
    input:
        tree = "results/tree.nwk",
        aln = "results/alignment.fasta",
        translations = expand("results/translation_{gene}.fasta", gene=genes),
        root = "config/reference.fasta",
        annotation = "config/reference.gb"
    output:
        node_data="results/muts.json"
    params:
        genes = genes,
        translations = "results/translation_%GENE.fasta"
    shell:
        """
        augur ancestral --tree {input.tree} --alignment {input.aln} \
                        --translations {params.translations} --genes {params.genes}\
                        --root-sequence {input.root} --annotation {input.annotation}\
                        --output-node-data {output.node_data}
        """

rule add_metadata:
    input:
        tree = "results/tree.nwk",
        metadata = "results/metadata.tsv"
    output:
        "results/node_metadata.json"
    run:
        import pandas as pd
        import json
        from Bio import Phylo   
        metadata = pd.read_csv(input.metadata, sep='\t', index_col='strain')
        metadata = metadata.to_dict(orient='index')
        T = Phylo.read(input.tree, "newick")
        nodes = {}
        for node in T.get_terminals():
            nodes[node.name] = {k1:metadata.get(node.name)[k2] for k1,k2 in 
                                zip(['country', 'date', 'subtype', 'Strain_name'], 
                                    ['country', 'date', 'subtype', 'name'])}
        with open(output[0], "w") as f:
            json.dump({"nodes":nodes}, f)


rule export:
    input:
        tree = "results/tree.nwk",
        node_data = ["results/branch_lengths.json", "results/clades.json",
                     "results/muts.json", "results/node_metadata.json"],
        auspice_config = "config/auspice_config.json",
        metadata = "results/metadata.tsv"
    output:
        "results/tree.json"
    shell:
        """
        augur export v2 --tree {input.tree} --node-data {input.node_data} --metadata {input.metadata} \
                        --auspice-config {input.auspice_config} --output {output} --minify-json \
                        --include-root-sequence-inline
        """

rule clean:
    shell:
        """
        rm -rf results
        rm -r dataset/tree.json
        """