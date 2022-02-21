# Fastp, um pré-processador FASTQ ultrarrápido com controle de qualidade
Autor: Luiz Carlos
Data: 10/01/2022

## Introdução

Geralmente temos duas abordagens diferentes ao aparar o adaptador:

1- usar uma ferramenta que pega um adaptador ou lista de adaptadores e os remove de cada sequência lida(trimmomatic, cutadapt).

2 - usar uma ferramenta que prevê adaptadores e os remove de cada sequência lida (fastp).

Ao lidar com dados públicos,a segunda abordagem pode ser empregada, ou seja, usar uma ferramenta que prevê adaptadores, 
pois possívelment não teremos as informações de preparo das bibliotecas.


## Fastp

* é um pré-processador FASTQ ultrarrápido com controle de qualidade útil e recursos de filtragem de dados. 
* realiza controle de qualidade, corte de adaptador, filtragem de qualidade, corte de qualidade por leitura 
e muitas outras operações com uma única varredura dos dados FASTQ. 
* Esta ferramenta foi desenvolvida em C ++ e possui suporte a multi-threading.

Opções:  

	-q, --qualified_quality_phred --> o valor de qualidade phred do qual uma base é qualificada. O padrão é 15, significando que apenas as bases com valor acima passam.   
	-l, --length_required ----> reads menores do que o tamanho especificado serão descartadas, o padrão é 15.  
	-f, --trim_front1 --------> Determina o número de bases a serem trimadas da extremidade 5" (início) da read1, default is 0;  
	-t, --trim_tail1 ---------> Determina o número de bases a serem trimadas da extremidade 3" (final) da read1, default is 0;  
	-F, --trim_front2 --------> Determina o número de bases a serem trimadas da extremidade 5" (início) da read2, default is 0; Se não especifidcado segue os valores da read1.  
	-T, --trim_tail2  --------> Determina o número de bases a serem trimadas da extremidade 3" (final) da read2, default is 0; Se não especifidcado segue os valores da read1.
	--detect_adapter_for_pe --> Especifica que estamos lidando com dados pareados.
	--overrepresentation_analysis -> Analisa a coleção de sequências para sequências que aparecem com muita frequência.
	--correction -------------> Tentará corrigir as bases com base em uma análise de sobreposição de read1 e read2.
	--cut_right --------------> Usará o corte de qualidade e varrerá a read do início ao fim em uma janela. Se a qualidade na janela estiver abaixo do necessário, a janela mais toda a sequência até o final é descartada e a read é mantida se ainda for longa o suficiente.  


## Executando fastp

Executando o fastp com as opções padrões: -q 15 e -l 15
Por padrão o fastp também procurar e reconhece os adptadores automaticamente.
```bash
fastp -i raw_fastq/reads_R1.fastq.gz -I raw_fastq/reads_R2.fastq.gz -o passedQC/reads_R1.fastq.gz -O passedQC/reads_R2.fastq.gz
```

Alterando valores de qualidade e comprimento das reads
```bash
fastp -q 30 -l 40 -i raw_fastq/reads_R1.fastq.gz -I raw_fastq/reads_R2.fastq.gz -o passedQC/reads_R1.fastq.gz -O passedQC/reads_R2.fastq.gz
```


Cortando um número fixo de nucleotídeos do início e final das reads
```bash
fastp -q 30 -l 40 -f 6 -F 6 -t 9 -T 9 -i raw_fastq/reads_R1.fastq.gz -I raw_fastq/reads_R2.fastq.gz -o passedQC/reads_R1.fastq.gz -O passedQC/reads_R2.fastq.gz
```


```bash
fastp --detect_adapter_for_pe
        --overrepresentation_analysis
        --correction --cut_right
        -i raw_fastq/reads_R1.fastq.gz -I raw_fastq/reads_R2.fastq.gz -o passedQC/reads_R1.fastq.gz -O passedQC/reads_R2.fastq.gz
```

## Referências:

https://github.com/OpenGene/fastp  
