IDS = ["Adl_Ss-w", "Dra_Ss-w", "Enn_Ss-w", "Gau_Ss-w"]
Assembler = ["norgal", "MitoZ", "MitoFlex", "GetOrganelle", "Novoplasty", "MITObim"]


rule all:
    input:
        expand("QUAST/report.tsv", id = IDS, assembler = Assembler)
        #expand("assemblies/norgal/Ecr_{id}/Ecr_{id}_norgal_genome.fa", id =IDS, assembler = Assembler)

#rule norgal_sub:
   # input:
     #   f = "raw_mt_reads/Ecr_{id}_1.fastq.gz",
     #   r = "raw_mt_reads/Ecr_{id}_2.fastq.gz"
    #output:
     #   f = "sub10mil/Ecr_{id}_1_sub.fastq.gz",
    #    r = "sub10mil/Ecr_{id}_2_sub.fastq.gz"
   # conda:
    #    "envs/seqtk.yml"
   # shell:
       # """
       # seqtk sample -s100 {input.f} 10000000 > {output.f} 
      #  seqtk sample -s100 {input.r} 10000000 > {output.r}
     #   """

#rule gunzip_sub:
    #input:
      #  f = "sub10mil/Ecr_{id}_1_sub.fastq.gz",
     #   r = "sub10mil/Ecr_{id}_2_sub.fastq.gz"
    #output:
     #   f = "sub10mil/Ecr_{id}_1_sub.fastq",
    #    r = "sub10mil/Ecr_{id}_2_sub.fastq"
   # shell:
    #    "gunzip {input.f} {input.r} > {output.f} {output.r}"


#rule norgal:
   # input:
      #  f = "sub10mil/Ecr_{id}_1_sub.fastq",
     #   r = "sub10mil/Ecr_{id}_2_sub.fastq"
    #output:
        #"assemblies/{assembler}/Ecr_{id}/Ecr_{id}_{assembler}_genome.fa"
    #shell:
        #"scripts/norgal.py -i {input.f} {input.r} -o {output} --blast"


rule trimmomatic:
    input:
        f = "raw_mt_reads/Ecr_{id}_1.fastq.gz",
        r = "raw_mt_reads/Ecr_{id}_2.fastq.gz"
    output:
        fout = "trimmed/{id}_1P_trim.fastq",
        funp = "trimmed/{id}_1P_unpaired.fastq",
        rout = "trimmed/{id}_2P_trim.fastq",
        runp = "trimmed/{id}_2P_unpaired.fastq"
    threads: 24
    conda:
        "envs/trimmomatic.yml"
    shell:
        "trimmomatic PE -threads {threads} {input.f} {input.r} {output.fout} {output.funp} {output.rout} {output.runp} ILLUMINACLIP:adapterseq/Adapters_PE.fa:2:30:10: LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:80"
  

rule download_GO_database:
    output:
       "get_organelle.db.ok"
    conda:
       "envs/getorganelle.yml" 
    shell:
       """
       get_organelle_config.py --clean
       get_organelle_config.py -a animal_mt
       touch {output}
       """


rule get_organelle:
    input:
        ok = rules.download_GO_database.output,
        f = "trimmed/{id}_1P_trim.fastq",
        r = "trimmed/{id}_2P_trim.fastq"
    output:
        "assemblies/{assembler}/Ecr_{id}/animal_mt.K115.complete.graph1.1.path_sequence.fasta"
    params:
        outdir = "assemblies/{assembler}/Ecr_{id}/"
    conda:
        "envs/getorganelle.yml"
    shell:
        "get_organelle_from_reads.py -1 {input.f} -2 {input.r} -o {params.outdir} -F animal_mt -R 10 -s seeds/Diphyllobothrium_stemmacephalum_mtgenome_NC_035881.1.fasta"

rule mitofelx:
    input:
        f = "trimmed/{id}_1P_trim.fastq",
        r = "trimmed/{id}_2P_trim.fastq"
    output:
        "assemblies/{assembler}/Ecr_{id}/Ecr_{id}.picked.fa"
    params:
        outdir = "assemblies/{assembler}/Ecr_{id}/"
    conda:
        "envs/mitoflex.yml"
    shell:
        "scripts/MitoFlex.py all --workname Ecr_{id} --threads 8 --insert-size 167 --fastq1 {input.f} --fastq2 {input.r} --genetic-code 14 --clade Platyhelminthes 1>m.log 2>m.err"

#rule NOVOplasty:
#    input:
#       "NOVOconfig_Ecr_{id}.txt"
#    output: 
#        "assemblies/{assembler}/Ecr_{id}/Ecr_{id}_novoplasty.fasta"
#    conda:
#        "envs/novoplasty.yml"
#    shell:
#        "scripts/NOVOPlasty4.2.1.pl -c {input}"

rule interleave:
    input:
        f = "trimmed/{id}_1P_trim.fastq",
        r = "trimmed/{id}_2P_trim.fastq"
    output:
        "trimmed/{id}_interleaved_trim.fastq"
    conda:
        "envs/bbmap.yml"
    shell:
        "reformat.sh in1={input.f} in2={input.r} out={output}"

rule MITObim:
    input:
        "trimmed/{id}_interleaved_trim.fastq"
    output:
        "assemblies/{assembler}/Ecr_{id}/Ecr_{id}_{assembler}.fasta"
    params:
        id = "{id}"
    singularity:
        "docker://chrishah/mitobim:v.1.9.1"
    shell:
        "scripts/MITObim.pl -sample Ecr_{params.id} -ref Diphyllobothrium_stemmacephalum_mtgenome -readpool trimmed/{params.id}_interleaved_trim.fastq --quick seeds/Diphyllobothrium_stemmacephalum_mtgenome_NC_035881.1.fasta -end 100 --denovo --paired --clean --NFS_warn_only &> log"
         

rule quast:
    input:
        #norgal = expand("assemblies/{assembler}/Ecr_{id}/Ecr_{id}_{assembler}_genome.fa", id = IDS,  assembler = Assembler[0]),
        #MitoZ = expand("assemblies/{id}_{assembler}/{id}_{assembler}-scaffolds.fa", id = IDS, assembler = Assembler),
        MitoFlex = expand("assemblies/{assembler}/Ecr_{id}/Ecr_{id}.picked.fa", id = IDS, assembler = Assembler[2]),
        GetOrganelle = expand("assemblies/{assembler}/Ecr_{id}/animal_mt.K115.complete.graph1.1.path_sequence.fasta", id = IDS, assembler = Assembler[3]),
        #Novoplasty = expand("assemblies/{assembler}/Ecr_{id}/Ecr_{id}_novoplasty.fasta", id = IDS, assembler = Assembler[4]),
        MITObim = expand("assemblies/{assembler}/Ecr_{id}/Ecr_{id}_{assembler}.fasta", id = IDS, assembler = Assembler[5])
    output:
        "QUAST/report.tsv"
    params:
        outdir = "/QUAST/",
    conda:
        "envs/quast.yml"
    shell:
        "quast -o {params.outdir} {input.GetOrganelle} {input.MitoFlex} {input.MITObim}"
