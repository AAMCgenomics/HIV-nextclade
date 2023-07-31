
genes = ["pro", "RT-p66", "RT-p51", "INT",
         "p17", "p24", "p2", "p7", "p1", "p6",
         "vif", "vpr", "tat",  "rev", "vpu",
         "gp120", "gp41", "nef"]

ruleorder: group_trees>tree

rule assemble:
    input:
        "config/reference.fasta",
        "config/annotation.gff",
        "config/tag.json",
        "config/qc.json",
        "config/primers.csv",
        "config/virus_properties.json",
        "results/tree.json"
    output:
        directory("dataset")
    shell:
        """
        mkdir -p {output}
        cp {input} {output}
        """

rule name_by_accession:
    input:
        "data/sequences.fasta",
        "data/additional.fasta"
    output:
        "data/sequences_renamed.fasta"
    run:
        from Bio import SeqIO
        accessions = set()
        with open(output[0], "w") as f:
            for fname in input:
                for record in SeqIO.parse(fname, "fasta"):
                    record.id = record.id.split(".")[-1]
                    if record.id in accessions:
                        continue
                    accessions.add(record.id)
                    record.description = ""
                    SeqIO.write(record, f, "fasta")

rule align:
    input:
        sequences = "data/sequences_renamed.fasta",
        ref = "config/reference.fasta",
        annotation = "config/annotation.gff"
    output:
        "results/alignment.fasta"
    params:
        translations = lambda w:"results/translation_{gene}.fasta"
    threads: 84
    shell:
        """
        ./nextalign run -j {threads} -r {input.ref} -m {input.annotation} --output-fasta {output} --output-translations {params.translations} \\
                        --output-insertions results/insertions.csv --include-reference --gap-alignment-side right \\
                        --penalty-gap-open 12 --penalty-gap-open-out-of-frame 14 --penalty-gap-open-in-frame 10 \\
                        {input.sequences}
        """

rule make_metadata:
    input:
        "data/sequences.fasta", "data/additional.fasta"
    output:
        "results/metadata.tsv"
    run:
        from Bio import SeqIO
        import pandas as pd
        metadata = [{'strain':'NC_001802','subtype':'B', 'country':'FR', 'date':1983, 'accession':'NC_001802', 'name':"HXB2_reference"}]
        fields = {'strain': -1, 'subtype': 1, 'country': 2, 'date': 3, 'accession': -1,  'name': [4,-1]}
        accessions = {metadata[0]['strain']}
        for fname in input:
            for record in SeqIO.parse(fname, "fasta"):
                entries = record.id.split(".")
                if len(entries)<6:
                    print(entries)
                    continue
                datum = {k: entries[v] if type(v)==int else '.'.join(entries[v[0]:v[1]]) for k, v in fields.items()}
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

rule split_by_subtype:
    input:
        alignment = "results/alignment.fasta",
        metadata = "results/metadata.tsv"
    output:
        alignments = directory("results/subtype_alignments")
    run:
        import pandas as pd
        import os
        from Bio import SeqIO
        os.makedirs(output.alignments)
        metadata = pd.read_csv(input.metadata, sep='\t', index_col='strain')
        for subtype, count in metadata.subtype.value_counts().items():
            if count<3: continue
            with open(output.alignments + f"/{subtype}.fasta", "w") as f:
                for record in SeqIO.parse(input.alignment, "fasta"):
                    if metadata.loc[record.id,'subtype'] == subtype:
                        SeqIO.write(record, f, "fasta")


rule trees_by_subtype:
    input:
        alignments = directory("results/subtype_alignments"),
    output:
        trees = directory("results/subtype_trees")
    shell:
        """
        mkdir -p {output.trees}
        for f in {input.alignments}/*.fasta; do
            augur tree --alignment $f --output {output.trees}/$(basename $f .fasta).nwk
        done
        """

rule group_trees:
    input:
        trees = directory("results/subtype_trees"),
    output:
        "results/tree_raw.nwk"
    run:
        from Bio import Phylo
        import glob


        T = Phylo.BaseTree.Tree()
        for t in glob.glob(input.trees + "/*.nwk"):
            sub_tree = Phylo.read(t, "newick")
            sub_tree.root_at_midpoint()
            sub_tree.root.branch_length = 0.1
            T.root.clades.append(sub_tree.root)
        Phylo.write(T, output[0], "newick")




rule splice_alignment_for_tree:
    input:
        "results/alignment.fasta"
    output:
        "results/alignment_spliced.fasta"
    run:
        from Bio import SeqIO
        with open(output[0], "w") as f:
            for record in SeqIO.parse(input[0], "fasta"):
                record.seq = record.seq[1798:3415]
                SeqIO.write(record, f, "fasta")

rule tree:
    input:
        aln = "results/alignment_spliced.fasta"
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree --alignment {input.aln} --output {output.tree}
        """

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
        annotation = "config/annotation.json"
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

rule export:
    input:
        tree = "results/tree.nwk",
        node_data = ["results/branch_lengths.json", "results/clades.json",
                     "results/muts.json"],
        auspice_config = "config/auspice_config.json"
    output:
        "results/tree.json"
    shell:
        """
        augur export v2 --tree {input.tree} --node-data {input.node_data} \
                        --auspice-config {input.auspice_config} --output {output} --minify-json
        """
