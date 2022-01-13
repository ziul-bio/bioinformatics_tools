# Mapeamento com HISAT2

HISAT é um programa de alinhamento rápido e sensível. Com ele foi desenvolvido um novo esquema de indexação baseado na transformação
de Burrows-Wheeler (BWT) e o índice FM, chamado de indexação hierárquica, que emprega dois tipos de índices, FM global, representando todo o genoma, e
índices FM locais, separados por pequenas regiões cobrindo coletivamente o genoma. 

O índice hierárquico para o genoma humano (cerca de 3 bilhões de pb) inclui ~ 48.000 índices FM locais, cada um representando uma região genômica de ~ 64.000 pb.
Como base para o alinhamento sem lacunas, o índice FM é extremamente rápido com um baixo uso de memória, conforme demonstrado por Bowtie. Além disso, o HISAT 
fornece várias estratégias de alinhamento projetadas especificamente para mapear diferentes tipos de leituras de RNA-seq. Juntos, o HISAT permite um alinhamento
extremamente rápido e sensível de leituras, em particular aquelas que abrangem dois exons ou mais. Como resultado, o HISAT é muito mais rápido e apresenta um boa
qualidade de alinhamento. 

Embora use um grande número de índices, o requisito de memória do HISAT ainda é modesto, aproximadamente 4,3 GB para o genoma humanos. 
HISAT usa a implementação Bowtie2 para lidar com a maioria das operações no índice FM. Além do alinhamento emendado, o HISAT lida com leituras envolvendo indels 
e suporta um modo de alinhamento de extremidades emparelhadas. Vários processadores podem ser usados simultaneamente para obter maior velocidade de alinhamento. 
O HISAT produz alinhamentos no formato SAM, permitindo a interoperação com um grande número de outras ferramentas que usam SAM. 
HISAT é distribuído sob a licença GPLv3 e é executado na linha de comando em Linux, Mac OS X e Windows.
[link](http://daehwankimlab.github.io/hisat2/)


### Install hisat2
```bash
conda install -c bioconda hisat2
```


#### Opções do HISAT2:
	-p # threads  
	-t # time  
	-f # arquivo fasta como arquivo de entrada  
	-x # reference genome  
	-x # reference genome  
	-1 # first mate  
	-2 # seconde mate  
	-S # sample output  
	--dta # Esta opção gera um relatório feito sob medida para os programas de reconstrução, incluindo StringTie.   
	--rna-strandness # Especifique as informações de orientação da fita: o padrão é unstranded.
	-- dta-cufflinks parameter makes hisat2 report alignment tailored specifically for Cufflinks.

A opção - downstream-transcriptome-assembly (--dta), é necessária caso se queira realizar "transcriptome assembly".
A opção requer comprimentos de âncora mais longos para a descoberta de novos transcritos, levando a menos alinhamentos com short-anchors, o que ajuda os 
"transcriptome assembly" a melhorar significativamente o uso de memória.

## Determinando a direcionalidade da biblioteca (- -rna-strandness)

Para que seja escolhido o parâmetro correto, é preciso descobrir qual foi método usado na preparar as bibliotecas. Para isso, é preciso testar a orientação/"strandness" das reads presentes nos arquivos fastq de RNA-Seq. 

Existem pelo menos dois métodos principais para a síntese de cDNA em bibliotecas "stranded":

* o método dUTP, que preserva a fita complementar

* o método Illumina, que liga diretamente ao RNA os ligantes, preserva a orientação original do RNA

É muito provável que o bioinformata não receba as informações da preparação da biblioteca de cDNA o que pode dificultar a escolha
dos parâmetros corretos. Por isso, realizar este teste agora, evita problemas futuros no decorrer da pipeline, como:

* A tag XS no arquivo BAM resultante do alinhamento, conterá informações incorretas da fita. A tag XS é usada por programas de montagem de transcrição, como Cufflinks e Stringtie, e também Cuffdiff a usa.

* Contagens incorretas de reads para alguns genes ou deixar de contar reads válidas para outros, na etapa de featureCounts.

[Mais informações sobre os tipos de bibliotecas](https://chipster.csc.fi/manual/library-type-summary.html)

***Obs: O pacote how_are_we_stranded_here.py em Python é usado para testar a orientação ou "strandness" das reads presentes nos arquivos fastq de RNA-Seq.
Porém, tive problemas na execução do script, mesmo usando a versão recomendado do kallisto, segue o link: [how_are_we_stranded_here](https://github.com/betsig/how_are_we_stranded_here)***


### Determinando --rna-strandness

Abaixo descrevo e executo os passos que o script how_are_we_stranded_here.py realizaria.

Requisitos:

* kallisto == 0,44.x

* RSeQC

* Arquivo de anotação de transcriptoma (arquivo .fasta - por exemplo, .cdna.fasta do ensembl)

* Arquivo de anotação no formato gtf, correspondente ao organismo em estudo.

* Arquivo gtf convertido em bed12. (Usei a ferramenta "gtf to bed12" do galaxy para fazer a converção. [link](https://usegalaxy.org/)) 

```bash
# criar um índice kallisto do transcriptoma do organismo em estudo.
kallisto index --index="Mus_musculus_index" reference_genome/Mus_musculus.GRCm39.cdna.all.fa

# mapear um subconjunto de leituras para o transcriptoma e usar o argumento --genomebam de kallisto, produzindo um arquivo pseudoa-lignment.bam ordenado.
kallisto quant --threads=4 --index=kallisto/Mus_musculus_index --bootstrap-samples=10 --output-dir kallisto/kallisto_results/ --genomebam --gtf reference_genome/Mus_musculus.GRCm39.104.gtf 132_S1/132_S1_L001_R1_001.fastq.gz 132_S1/132_S1_L001_R2_001.fastq.gz 

# Executar infer_experiment.py do RSeQC para verificar qual direção as reads estão alinhadas em relação à fita transcrita.
infer_experiment.py -r reference_genome/Mus_musculus.GRCm39.104.gtf.bed12 -i kallisto/kallisto_results/pseudoalignments.bam
```

Output:

	Reading reference gene model reference_genome/Mus_musculus.GRCm39.104.gtf.bed12 ... Done  
	Loading SAM/BAM file ...  Total 200000 usable reads were sampled  
	This is PairEnd Data  
	Fraction of reads failed to determine: 0.0095  
	Fraction of reads explained by "1++,1--,2+-,2-+": 0.0019  
	Fraction of reads explained by "1+-,1-+,2++,2--": 0.9886  

Interpretação:  

0.95% foram mapeadas em regiões do genoma que não foi possível determinar “standness of transcripts” (como regiões que tiveram ambas as fitas transcritas).  
Os outros 98.86% das reads, a grande maioria é explicada por, "1+-,1-+,2++,2--", sugerindo um dataset "strand-specific -reverse".

[Para mais informações sobre a interpretação na secção RSeQC](#rseqc-an-rnaseq-quality-control-package)


## Voltando para o HISAT2 - Contruindo um indíce genômico dos sitios de splicing e informações de exon

***Primeiro, usando os scripts python incluídos no pacote HISAT2, extraia o sitios de splicing e as informações do exon 
do arquivo de anotação de gene (GTF, mas caso queira, conseguirá fazer sem ele).***

```bash
extract_splice_sites.py genome.gtf > genome.ss
extract_exons.py genome.gtf > genome.exon
```

Opções:  

	--ss   path/para os arquivos .ss   
Obs: Esta opção deve ser usada com a opção --exon. Forneça a lista de splice sites (no formato HISAT2, 4 colunas)  

	--exon path/para os arquivos .exon  
Esta opção deve ser usada com a opção --ss acima. Forneça a lista de exons (no formato HISAT2, 3 colunas).  


## Criando um índice HISAT2 das regiões seleciondas:

Opinião: "Não consegui executar esta etapa, recebi a mensagem de erro: 'ran out of memory' "
```bash
### Opções:
--ss   path/para os arquivos .ss    # Esta opção deve ser usada com a opção --exon. Forneça a lista de splice sites (no formato HISAT2, 4 colunas)

--exon path/para os arquivos .exon  # Esta opção deve ser usada com a opção --ss acima. Forneça a lista de exons (no formato HISAT2, 3 colunas).  

hisat2-build --ss genome.ss --exon genome.exon genome.fasta genoma_index   
```


## Criando um índice do genoma de referência:

O genonas de referência para a espécie em estudo (usado para contruir o índice) pode ser encontrado no site ensembl [link](https://www.ensembl.org/index.html).
Ensembl é um navegador de genomas de vertebrados que oferece suporte à pesquisa em genômica comparativa, evolução, variação de sequência 
e regulação da transcrição. O Ensembl anota genes, calcula alinhamentos múltiplos, prediz a função regulatória e coleta dados de doenças.
As ferramentas de ensembl incluem BLAST, BLAT, BioMart e o Variant Effect Predictor (VEP) para todas as espécies com suporte.

```bash
hisat2-build -p 4 reference_genome.fa genome
```
O comando acima gera 8 arquivos .ht2, sendo eles usados no alinhamento das reads.


# Alinhamento com o genoma de referência

Os comando abaixo usam as informações que adquirimos anteriormente, como direcionalidade das bibliotecas e o índice genômico criado.
```bash
# Exemple 1, para paired-end reads, arquivo de saída sam:
hisat2 -p 5 -t --dta --rna-strandness RF -x .../reference_genome/genome -1 sample_R1.fq -2 sample_R2.fq -S sample.sam

# Exemple 2, para paired-end reads, direcionando o output para um arquivo bam ordenado:
hisat2 -p 5 -t --dta --rna-strandness RF -x .../reference_genome/genome  -1 .../1-Raw_data/fastq/SRR00000001_1.fastq -2 .../1-Raw_data/fastq/SRR00000001_2.fastq | samtools view -b -o S1_L001_sorted.bam - 

```
Output do mapeamento com hisat2: 
![output do mapeamento com hisat2](images/output_alignment.png)

[Informações sobre a interpretação dos resultados](https://www.biostars.org/p/395017/)  
[Mais informações sobre hisat](http://daehwankimlab.github.io/hisat2/manual/)



## Alinhamento com spliced (Spliced alignment)

	--dta/--downstream-transcriptome-assembly = Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

	--dta-cufflinks = Report alignments tailored specifically for Cufflinks. In addition to what HISAT2 does with the above option (–dta), With this option, HISAT2 looks for novel splice sites with three signals (GT/AG, GC/AG, AT/AC), but all user-provided splice sites are used irrespective of their signals. HISAT2 produces an optional field, XS:A:[+-], for every spliced alignment.

	--novel-splicesite-outfile = In this mode, HISAT2 reports a list of splice sites in the file:
	chromosome name tab genomic position of the flanking base on the left side of an intron tab genomic position of the flanking base on the right tab strand (+, -, and .)
	’.’ indicates an unknown strand for non-canonical splice sites.


```bash
# Para criar um arquivo specífico para uso com cufflink:
hisat2 --dta-cufflinks -p 5 -t -x ../reference/genome.fa -1 ../sample_R1_1.fq -2 ../sampl_R1_2.fq -S sample_R1.sam --rna-strandness RF 

hisat2 --dta-cufflinks -p 5 -t -x ../reference/genome.fa -1 ../sample_R2_1.fq -2 ../sample_R2_2.fq -S sample_R2.sam --rna-strandness RF


# Para criar um arquivo splice junction:
hisat2 --dta -p 5 -t -x ../reference/genome.fa -1 ../sample_R1_1.fq -2 ../sampl_R1_2.fq -S sample_R1.sam --novel-splicesite-outfile sample_R1.junctions --rna-strandness RF 

hisat2 --dta -p 5 -t -x ../reference/genome.fa -1 ../sample_R2_1.fq -2 ../sample_R2_2.fq -S sample_R2.sam --novel-splicesite-outfile sample_R2.junctions --rna-strandness RF
```


# Alternativa de alinhamento em caso de fastq com multiplas lanes

##### Primeira opção:
Realiza o alinhamento de cada lane, e ao final se realiza o merge dos arquivos bam:
```bash
## Uso do pipe para direcionar o output do alinhamento para um arquivo BAM, sem a necessidade do intermediário SAM. Estes arquivos ocupam muito espaço.

# Nesta caso será necessário ordenar cada BAM após a conversão. 
hisat2 --dta -p 5 -t -x ~/bioinfo/rna_seq/1-raw_data/reference_genome/genome -1 S1_L001_R1.fastq -2 S1_L001_R2.fastq | samtools view -b -o S1_L001.bam -

# Direcionando o output do alinhamento direto para um arquivo bam ordenado. Usando a função samtools sort: "O hífen(-) ao final do comando representa 
# informa ao comando para usar como input presente no pipe, neste caso o output do alinhamento"
hisat2 --dta -p 5 -t -x ~/bioinfo/rna_seq/1-raw_data/reference_genome/genome -1 S1_L001_R1.fastq -2 S1_L001_R2.fastq | samtools sort -o S1_L001.bam -

# Direcionando o output para bam ordenado, mas mantendo o cabeçalho (header) bam.
hisat2 --dta -p 5 -t -x ~/bioinfo/rna_seq/1-raw_data/reference_genome/genome -1 S1_L001_R1.fastq -2 S1_L001_R2.fastq | samtools view -h -b - | samtools sort -o S1_L001.bam
```

Unindo os arquivos bam, resultantes de cada alinhamento
```bash
# Merge, usado para arquivos sam, bam e cram, toma como entrada arquivos ordenados e gera um arquivo tbm ordenado. *** estou usando esta opção.***
samtools merge S1_merged.bam S1_L001.bam S1_L002.bam S1_L003.bam S1_L004.bam

# cat, usado para arquivos bam e cram, e o dicionário de sequência dos arquivos sendo concatenados precisa ser idêntico
samtools cat 
``` 
Média 40min por lane. Sendo o merge mais uns 30min.


##### Segunda opção:

Passar uma lista separada por vírgulas ao hisat2, com os arquivos fastqc(L001_R1.fastqc, L002_R1 ...) para o mate1 e uma outra lista, com os arquivos fastqc(L001_R2.fastqc, L002_R2 ...) para o mate2. Sendo o sort e merge realizado em pipe:
```bash
# Passar uma lista dos arquivo correspondente a cada lane separado por vírgula, e direcionar o output para um arquivo convertido à BAM ordenado, 
# OBS: O arquivo de output será muito grande, por isso já realizei a compactação. 
hisat2 -p 5 -t --dta --rna-strandness RF -x ~/bioinfo/rna_seq/1-raw_data/Reference_genome/genome -1 S1_L001_R1.fastq,S1_L002_R1.fastq,S1_L003_R1.fastq,S1_L004_R1.fastq -2 S1_L001_R2.fastq,S1_L002_R2.fastq,S1_L003_R2.fastq,S1_L004_R2.fastq | samtools sort -o /mnt/d/HD_exter/S1_L001.bam -

hisat2 -p 5 -t --dta --rna-strandness RF -x ~/bioinfo/rna_seq/1-raw_data/reference_genome/genome -1 133_S2_L001_R1_001.fastq.gz,133_S2_L002_R1_001.fastq.gz,133_S2_L003_R1_001.fastq.gz,133_S2_L004_R1_001.fastq.gz -2 133_S2_L001_R2_001.fastq.gz,133_S2_L002_R2_001.fastq.gz,133_S2_L003_R2_001.fastq.gz,133_S2_L004_R2_001.fastq.gz | samtools sort -o ~/bioinfo/rna_seq/3-alignment_hisat2/133_S2.bam -

# Estou usando está opção, pois o sort junto com convert´precisa de uma versão do samtools igual ou superior 1.2. Além disto, a opção -h garante a inserção do head.
hisat2 -p 5 -t --dta --rna-strandness RF -x /mnt/c/Users/luiz_/wsl/rna_seq/1-raw_data/reference_genome/genome -1 133_S2_L001_R1_001.fastq.gz,133_S2_L002_R1_001.fastq.gz,133_S2_L003_R1_001.fastq.gz,133_S2_L004_R1_001.fastq.gz -2 133_S2_L001_R2_001.fastq.gz,133_S2_L002_R2_001.fastq.gz,133_S2_L003_R2_001.fastq.gz,133_S2_L004_R2_001.fastq.gz | samtools view -h -b - | samtools sort -o /mnt/c/Users/luiz_/wsl/rna_seq/3-alignment_hisat2/133_S2.bam
```
Média de 3:30min por amostra.

[Informações sobre multiplas lanes](#material-complementar)


## Referências:

https://anaconda.org/bioconda/hisat2  
http://www.htslib.org/doc/samtools.html   
http://rseqc.sourceforge.net/#infer-experiment-py   