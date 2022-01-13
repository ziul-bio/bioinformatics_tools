# Fastp, um pré-processador FASTQ ultrarrápido com controle de qualidade
Autor: Luiz Carlos
Data: 10/01/2022

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


## Referências:

https://github.com/OpenGene/fastp  
