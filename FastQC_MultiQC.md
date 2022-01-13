# Controle de qualidade de arquivos fastq com fastrqc e Multiqc

Os dados gerados a partir de tecnologias de sequenciamento podem gerar um grande número de leituras de 
sequência em um único experimento. No entanto, cada instrumento irá gerar diferentes tipos de erros, 
como sequência incorreta de nucleotídeos. Essas bases erroneamente chamadas se devem às limitações técnicas 
de cada plataforma de sequenciamento.

Portanto, é necessário identificar e excluir tipos de erros que podem impactar a interpretação da análise a jusante.

O controle de qualidade da sequência é, portanto, um primeiro passo essencial em sua análise, evitando dores de cabeça posteriormente.


## Fastqc

FastQC é um programa de controle de qualidade que nos permite realizar várias verificações de controle de qualidade
em dados brutos de sequênciamento gerados por NGS, como as plataformas Illumina e ABI SOLiD no formato FASTQ. 
Gera como arquivo de saída um relatório abrangente de várias páginas sobre a composição e qualidade das leituras em formato HTML.


### Opções do fastqc:

    -o  = Diretório externo, que se deseja mover os reports.  
    -t  = Número de threads (núcleos de cpu)  
    -q  = "quiet" suprime todo o progresso, apenas reportando erros.  

Criando um novo diretório para enviar os reports:
```bash
mkdir ~/bioinfo/rna_seq/fastqc
```

Executando o fastqc e gerando os reports: 
```bash
fastqc -t 5 *.fastq # o símbolo * indica todos os arquivos que terminem em .fastq
```

Executando o fastqc e movendo os reports para o directório especificado
```bash
fastqc -t 5 *.fastq -o ~/bioinfo/rna_seq/2-CQ/fastqc
```

Movendo os reports após a execução do fastqc:
```bash
mv *fastqc* ~/bioinfo/rna_seq/2-QC/fastqc
```


## Multiqc

MultiQC é uma ferramenta de relatório que analisa estatísticas resumidas e arquivos de log gerados por outras ferramentas de bioinformática. 
O MultiQC foi projetado para ser colocado no final de pipelines de análises ou para ser executado manualmente ao termino da execução de algum programa.

Ao executar o MultiQC, ele pesquisa recursivamente em qualquer caminho de arquivo fornecido e ao encontra os arquivos que reconhece, analisa as informações relevantes
e gera um único arquivo de relatório HTML independente, além de criar um diretório "multiqc_data" onde armazena os dados analisados, para uso posterior.

### Opções:

    .            = O ponto serve para executar o multiqc no diretório de trabalho.  
    -f           = Sobrescreve reports anteriores, em a necessidade de --force.  
    --force      = Força e execução do programa.  
    -x, --ignore = ignorar arquivos no diretório.  
    -p, --export = Para salvar os plots em arquivos separados no report, no diretório /multiqc_plots como .png, .svg e .pdf.  
    -o, --outdir = Gera todos os reports no directório externo especificado.  

Obs: O multiqc usa os reports do fastqc como arquivos de entrada, por isso precisa ser executado no diretório onde tais arquivos estão.
Além disso, o MultiQC foi projetado para ser executados em diretórios, ao invés de arquivos. Por isso, é mais fácil e recomendável informar 
apenas o/os diretório(s) a serem analisados pelo programa.  


### Executando o multiqc
```bash
# Executando multiqc no diretório de trabalho.
multiqc .

# Especificando um diretório.
multiqc ~/bioinfo/rna_seq/2-QC/fastqc

# Especificando multiplos diretórios
multiqc ~/bioinfo/rna_seq/2-QC/ # esta opção o multiqc procurará por arquivos em todas as pastas.
multiqc ~/bioinfo/rna_seq/2-QC/fastqc/ -o ~/bioinfo/rna_seq/2-QC/  

# Especificando diretório e arquivos a serem analisados
multiqc data/*fastqc.gz

# Para forçar o programa Executar novamente.
multiqc --force . 

# Para ignorar algum arquivo
multiqc --ignore 'RSS00000001' .

# Gerando um relatório final
multiqc ~/bioinfo/rna_seq/2-QC/* ~/bioinfo/rna_seq/3-hitsat/* ~/bioinfo/rna_seq/4-featureCounts/* 

```
 
### Visualizando os reports
```bash
# Listando o diretório multiqc_data/
ls multiqc_data/

# Estou usando wsl-2 - Ubuntu, por isso criei um alias do wslwiew = open para abrir os arquivos html.
open multiqc_report.html
```


## Referências:

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
https://multiqc.info/docs/#running-multiqc