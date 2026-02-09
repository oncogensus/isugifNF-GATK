# isugifNF-GATK 
## (A descrição, instruções e adaptações do manual original de uso deste pipeline ainda estão em produção)

Para fazer a chamada de variantes de SVN, é utilizado o pipeline em [Nextflow](https://www.nextflow.io/) isugifNF/GATK. Este pip[eline automatiza o uso do Genome Analysis Toolkit (GATK), modificado para seguir a melhores práticas descritas no Bioinformatic Workbook: GATK Best Practices Workflow for DNA-Seq. Este pipeline estende a funcionalidade do GATK, que originalmente é focado em dados de DNA-Seq, para permitir a detecção de variantes também a partir de dados de RNA-Seq e de sequenciamento de leituras longas. Ele fornece uma solução abrangente para diferentes plataformas de sequenciamento, ampliando as capacidades de análise de variantes.

## Instalação

Você precisará de uma versão funcional do Nextflow. Veja [aqui](https://www.nextflow.io/docs/latest/getstarted.html#requirements)￼como instalar o Nextflow. Módulos do Nextflow estão disponíveis em alguns ambientes de computação HPC. Todo processo de download das ferramentas é automatizado pelo pipeline via Nextflow e Singularity.

<details><summary>Modulos de clusters HPC</summary>

```
# === Nova
module load gcc/7.3.0-xegsmw4 nextflow
module load singularity
NEXTFLOW=nextflow

# === Condo
module load gcc/7.3.0-xegsmw4 nextflow
module load singularity
NEXTFLOW=nextflow

# === Ceres
module load nextflow
# singularity already available, no need for module
NEXTFLOW=nextflow

# === Atlas (will need a local install of nextflow and will need the --account "projectname" flag)
module load singularity
NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow
```

</details>

Uma vez que o Nextflow esteja instalado, o comando abaixo pode ser usado para ver o menu do isugifNF/GATK

```
nextflow run isugifNF/GATK --help
```

<details><summary>See help statement</summary>

```
N E X T F L O W  ~  version 20.10.0
Launching `isugifNF/GATK` [big_kare] - revision: ca139b5b5f
   Usage:
   The typical command for running the pipeline are as follows:

   DNAseq:
     nextflow run main.nf --genome GENOME.fasta --reads "*_{R1,R2}.fastq.gz" --seq "dna" -profile singularity
     nextflow run main.nf --genome GENOME.fasta --reads_file READ_PATHS.txt --seq "dna" -profile singularity

   RNAseq:
     nextflow run main.nf --genome GENOME.fasta --gtf "genes.gtf" --reads "*_{R1,R2}.fastq.gz" --seq "rna" -profile singularity
  
   PacBio Long Reads:
     nextflow run main.nf --genome GENOME.fasta --long_reads "*.fastq.gz" --seq "longread" -profile singularity

   Mandatory arguments:
    --seq                   Specify input sequence type as 'dna', 'rna', or 'longread' [default:'dna'].
    --genome                Reference genome fasta file, against which reads will be mapped to find Variant sites

   Read input arguments:
    --reads                 Paired-end reads in fastq.gz format, will need to specify glob (e.g. "*_{R1,R2}.fastq.gz")
    --reads_file            Text file (tab delimited) with three columns [readname left_fastq.gz right_fastq.gz]. Will need full path for files.
    --long_reads            Long read file in fastq.gz format, will need to specify glob (e.g. "*.fastq.gz")

   Optional analysis arguments:
    --invariant             Output invariant sites [default:false]
    --gtf                   Gene Transfer Format file, only required for RNAseq input [default:false]
    --window                Window size passed to bedtools for parallel GATK Haplotype calls [default:100000]

   Optional configuration arguments:
    -profile                Configuration profile to use. Can use multiple (comma separated)
                            Available: local, slurm, singularity, docker [default:local]
    --container_img         Container image used for singularity and docker [default:'docker://ghcr.io/aseetharam/gatk:master']
    
   GATK:
    --gatk_app              Link to gatk executable [default: 'gatk']
    --java_options          Java options for gatk [default:'-Xmx80g -XX:+UseParallelGC -Djava.io.tmpdir=$TMPDIR']
    --gatk_cluster_options  GATK cluster options [default:'false']
    
   Aligners:
    --bwamem2_app           Link to bwamem2 executable [default: 'bwa-mem2']
    --star_app              Link to star executable [default: 'STAR']
    --star_index_params     Parameters for star index [default: ' ']
    --star_index_params     Parameters to pass to STAR index [default:'']
    --star_index_file       Optional: speedup by providing a prebuilt STAR indexed genome [default: 'false']
    --pbmm2_app             Link to pbmm2 executable [default: 'pbmm2']

   Other:
    --samtools_app          Link to samtools executable [default: 'samtools']
    --bedtools_app          Link to bedtools executable [default: 'bedtools']
    --datamash_app          Link to datamash executable [default: 'datamash']
    --vcftools_app          Link to vcftools executable [default: 'vcftools']

   Optional other arguments:
    --outdir                Output directory [default:'./GATK_Results']
    --threads               Threads per process [default:4 for local, 16 for slurm] 
    --queueSize             Maximum jobs to submit to slurm [default:40]
    --account               HPC account name for slurm sbatch, atlas and ceres requires this
    --help                  Print this help message

NOTE: Genome name should not contain any periods before the ".fasta"
```

</details>

### Container Singularity

As ferramentas necessárias para o workflow estão incluídas no container [aseetharam/gatk:latest](https://github.com/aseetharam/gatk)￼ e devem ser baixadas automaticamente pelo Nextflow. (Você só precisará executar singularity pull se a conexão com o site estiver instável.)

#### Para baixar a imagem do container Singularity

```
singularity pull --name gatk.sif shub://aseetharam/gatk:latest
```

#### No Nextflow use a flag `-with-singularity` para usar a imagem baixada.

```
nextflow run isugifNF/GATK \
  --genome "test-data/ref/b73_chr1_150000001-151000000.fasta" \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  -profile slurm \
  -with-singularity gatk.sif
```

</details>

### Dataset de teste

Um conjunto de dados simples para teste da ferramenta está disponível no [ISU Box](https://iastate.app.box.com/v/gatk-test-data)￼. Esse conjunto contém um genoma pequeno (uma porção do cromossomo 1, B73v5) e leituras curtas de Illumina para 26 linhagens NAM (incluindo B73) e para a linhagem B73Ab10 (total de 27 linhagens).

Apenas as leituras que mapeiam para a região do genoma v5 estão incluídas, para que o teste possa ser executado rapidamente.

Há exemplos com múltiplos arquivos pertencentes à mesma linhagem NAM, assim como exemplos com um único arquivo por linhagem NAM, para garantir que ambas as condições funcionem corretamente.

O arquivo VCF final deve conter exatamente 27 indivíduos (linhagens).

Para baixar o dataset de teste, use os seguintes comandos:
```
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz

ls -1 test-data/
#> fastq          # <= Pasta com a reads (paired)
#> read-group.txt
#> ref            # <= Para com o genoma de referência
```

## Running the Pipeline

This pipeline can process DNAseq, RNAseq and PacBio long read data, which require different arguments and different input files. Please jump to the section based off of your input below.

### GATK variant calling for DNAseq data

Because this pipeline was initially developed to process DNAseq data, we have provided a test DNAseq dataset. You can fetch the pipeline and the test-data folder.

```
# Fetch pipeline
nextflow run isugifNF/GATK --help

# Fetch the test-data folder from ISU box
wget https://iastate.box.com/shared/static/wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
tar -xf wt85l6s4nw4kycm2bo0gpgjq752osatu.gz
```

The general format of a run with the pipeline is to provide a genome file (`--genome`) and Illumina Paired-End Reads files (`--reads` or `--reads_file`).

```
nextflow run isugifNF/GATK \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads "test-data/fastq/*_{R1,R2}.fastq.gz" \
  -profile slurm,singularity \
  -resume
```

or

```
nextflow run isugifNF/GATK \
  --genome test-data/ref/b73_chr1_150000001-151000000.fasta \
  --reads_file read-path.txt \
  -profile slurm,singularity \
  -resume
```

If you are on a HPC (Nova/Condo/Ceres/Atlas), it is highly recommend to use the `submit_nf.slurm` script. The `# === Load Modules` section will need to be modified to get nextflow and singulariy running.

<details><summary>See Module Changes</summary>

```
# === Nova
module load gcc/7.3.0-xegsmw4 nextflow
module load singularity
NEXTFLOW=nextflow

# === Condo
module load gcc/7.3.0-xegsmw4 nextflow
module load singularity
NEXTFLOW=nextflow

# === Ceres
module load nextflow
# singularity already available, no need for module
NEXTFLOW=nextflow

# === Atlas
module load singularity
NEXTFLOW=/project/isu_gif_vrsc/programs/nextflow
```

</details>

### Using the --reads_file

Instead of using a pattern to specify paired reads (`"test-data/fastq/*_{R1,R2}.fastq.gz"`), we can use a tab-delimited textfile to specify the path to left and right files. The textfile should contain three columns: readname, left read path, right read path.

In our case, for the `test-data` run the following one-liner to generate a `read-path.txt`.

```
for f in test-data/fastq/*_R1.fastq.gz; do echo -e "$(basename $f |cut -f 1 -d "_")\t$(realpath $f)\t$(realpath $f | sed 's/_R1.fastq.gz/_R2.fastq.gz/g')"; done > read-path.txt
```

<details><summary>See example <b>read-path.txt</b></summary>

```
BioSample01	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample01_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample01_R2.fastq.gz
BioSample02	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample02_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample02_R2.fastq.gz
BioSample03	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample03_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample03_R2.fastq.gz
BioSample04	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample04_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample04_R2.fastq.gz
BioSample05	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample05_R1.fastq.gz	/Users/jenchang/Maize_WGS_Build/test-data/fastq/BioSample05_R2.fastq.gz
```

</details>

#### DNAseq final output

The final output will be in a `results` folder. SNPs will be in the VCF file, probably the file with the longest name (e.g. `first-round_merged_snps-only_snp-only.pass-only.vcf`).

```
results/
  |_ 01_MarkAdapters/            #<= folders contain intermediate files
  |_ 02_MapReads/
  |_ 03_PrepGATK/
  |_ 04_GATK/
  |_ 05_FilterSNPs/
  |  |_ first-round_merged_snps-only_sorted_snp-only.pass-only.vcf     #<= final SNP file
  |
  |_ report.html
  |_ timeline.html               # <= runtime information for all processes

```
## O arquivo de configuração do Nextflow nf_local.config
Esta configuração otimiza o pipeline isugifNF/GATK para execução em ambientes onde o armazenamento compartilhado (NFS) pode ser instável ou lento durante operações intensivas do Picard e GATK. Abaixo seguem algumas observações e explicações sobre o arquivo nf_local.config.

#### 1. Observabilidade e Rastreamento
As seções iniciais habilitam a geração de relatórios detalhados após a execução:

Report/Timeline: Gera arquivos HTML com o tempo de execução e uso de recursos por processo.

Trace: Fornece uma tabela detalhada (status, hash, memória real usada) de cada task.

#### 2. Gestão de Arquivos Temporários (Local TMP)
Esta é a parte mais crítica do script. Processos baseados em Java (Picard/GATK) utilizam o SortingCollection, que despeja arquivos imensos no diretório temporário quando a memória RAM não é suficiente para ordenar os dados.

Estratégia: Foi definida a variável LOCAL_TMP_BASE apontando para o /tmp local do nó de computação (disco rígido físico da máquina, não a rede).

Objetivo: Evitar erros de "Disk quota exceeded" ou "I/O Error" comuns quando o NFS tenta lidar com milhares de pequenos arquivos temporários.

#### 3. Configurações de Processo (Defaults vs. Específicos)
O arquivo divide os processos em duas categorias principais:

A. I/O-Bound (Limitados por Disco)
Processos como FastqToSam, MergeBamAlignment, MarkDuplicates e SortAndFixTags.

maxForks = 1: Força a execução de apenas uma amostra por vez para esses processos. Isso evita que múltiplas tarefas disputem a largura de banda do disco simultaneamente, o que causaria lentidão geral.

beforeScript: Cria e limpa manualmente pastas no /tmp local antes de iniciar o processo, garantindo que o Singularity/Apptainer e o Java usem esse espaço rápido.

Gestão de Heap Java (-Xmx): O script define o uso de memória Java de forma conservadora (geralmente menor que a memória total do container) para deixar margem para o sistema operacional e evitar o encerramento do processo pelo gerenciador de memória (OOM Killer).

B. CPU-Bound (Limitados por Processamento)
Processos como bwamem2 e HaplotypeCaller.

maxForks = 2: Permite um pouco mais de paralelismo, já que o gargalo aqui é o cálculo matemático e não apenas a escrita em disco.

Recursos Elevados: Alocação de até 16 CPUs e 96GB de RAM para acelerar o alinhamento e a chamada de variantes.

#### 4. Controle do Executor
queueSize = 50: Limita o Nextflow a submeter no máximo 50 tarefas totais para a fila do cluster ao mesmo tempo, mantendo o controle sobre a carga total no servidor.

```
report   { enabled = true; overwrite = true }
timeline { enabled = true; overwrite = true }
trace    { enabled = true }

// =====================================================
// TMP local (somente para processos Picard “sensíveis”)
// - Evita NFS flakey no SortingCollection.*.tmp
// - TMP fixo por processo (sem task.hash)
// =====================================================
def LOCAL_TMP_BASE = "/tmp/${System.getenv('USER') ?: 'user'}/isugif_tmp"

// -----------------------------------------------------
// Defaults (baixos) — o pipeline sobe nos processos pesados
// -----------------------------------------------------
process {
  cpus   = 2
  memory = 16.GB
  time   = '48h'

  // Ajuda a não entupir o storage com work dirs antigos
  scratch = false
}

// -----------------------------------------------------
// Java default conservador
// (os processos I/O-bound vão sobrescrever)
// -----------------------------------------------------
params {
  java_options = "-Xmx16g -XX:+UseParallelGC -Djava.io.tmpdir=${System.getenv('TMPDIR')}"
}

// ===============================
// I/O-bound (serializar mais)
// ===============================
//
// Aqui é onde seu erro acontece.
// Com muitos genomas, o gargalo é disco/temporários, não CPU/RAM.
//
process {

  // 1) PRINCIPAL CULPADO: FastqToSam
  withName: 'DNA_VARIANT_CALLING:FastqToSam' {
    cpus     = 2
    memory   = 48.GB
    time     = '48h'
    maxForks = 1

    // reduzir heap -> reduz spill / temp e volume de escrita
    ext.java_opts = "-Xmx16g -XX:+UseParallelGC -Djava.io.tmpdir=${System.getenv('TMPDIR')}"
  }

  // 2) MergeBamAlignment: usa SortingCollection -> precisa TMP local
  withName: 'DNA_VARIANT_CALLING:MergeBamAlignment_DNA' {
    cpus     = 3
    memory   = 64.GB
    time     = '72h'
    maxForks = 1

    // FORÇA TMP LOCAL (não-NFS) dentro do task
    beforeScript = """
      export TMPDIR="${LOCAL_TMP_BASE}/MergeBamAlignment_DNA"
      export SINGULARITY_TMPDIR="\$TMPDIR"
      export APPTAINER_TMPDIR="\$TMPDIR"
      mkdir -p "\$TMPDIR" && chmod 700 "\$TMPDIR"
    """

    // garante java.io.tmpdir = TMPDIR local
    ext.java_opts = '-Xmx24g -XX:+UseParallelGC -Djava.io.tmpdir=$TMPDIR'
  }

  // SortAndFixTags também pode usar temp pesado (Picard-ish)
  withName: 'DNA_VARIANT_CALLING:SortAndFixTags_DNA' {
    cpus     = 3
    memory   = 64.GB
    time     = '72h'
    maxForks = 1

    beforeScript = """
      export TMPDIR="${LOCAL_TMP_BASE}/SortAndFixTags_DNA"
      export SINGULARITY_TMPDIR="\$TMPDIR"
      export APPTAINER_TMPDIR="\$TMPDIR"
      mkdir -p "\$TMPDIR" && chmod 700 "\$TMPDIR"
    """

    ext.java_opts = '-Xmx24g -XX:+UseParallelGC -Djava.io.tmpdir=$TMPDIR'
  }

  // MarkDuplicates costuma gerar MUITO temp e I/O (SortingCollection)
  withName: 'DNA_VARIANT_CALLING:MarkDuplicates_DNA' {
    cpus     = 4
    memory   = 96.GB
    time     = '96h'
    maxForks = 1

    beforeScript = """
      export TMPDIR="${LOCAL_TMP_BASE}/MarkDuplicates_DNA"
      export SINGULARITY_TMPDIR="\$TMPDIR"
      export APPTAINER_TMPDIR="\$TMPDIR"
      mkdir -p "\$TMPDIR" && chmod 700 "\$TMPDIR"
    """

    ext.java_opts = '-Xmx32g -XX:+UseParallelGC -Djava.io.tmpdir=$TMPDIR'
  }

  // Se existir e for Java/GATK, também é I/O pesado
  withName: 'DNA_VARIANT_CALLING:SamToFastq_DNA' {
    cpus     = 2
    memory   = 32.GB
    time     = '48h'
    maxForks = 1
  }

  withName: 'DNA_VARIANT_CALLING:MarkIlluminaAdapters' {
    cpus     = 2
    memory   = 32.GB
    time     = '48h'
    maxForks = 1
  }
}

// ===============================
// CPU-bound (limitar para não “derreter” o disco)
// ===============================
process {

  withName: 'DNA_VARIANT_CALLING:bwamem2_mem' {
    cpus     = 16
    memory   = 96.GB
    time     = '72h'
    maxForks = 2
  }

  withName: 'DNA_VARIANT_CALLING:gatk_HaplotypeCaller_DNA' {
    cpus     = 12
    memory   = 96.GB
    time     = '96h'
    maxForks = 2
    ext.java_opts = "-Xmx48g -XX:+UseParallelGC -Djava.io.tmpdir=${System.getenv('TMPDIR')}"
  }
}

// ===============================
// Segurança extra: limite global de tasks concorrentes
// ===============================
executor {
  queueSize = 50
}
```
## Guia de Manutenção e Limpeza de Storage - Oncogensus

Este documento descreve as políticas de retenção de dados e limpeza para execuções do pipeline `isugifNF/GATK` no cluster. O objetivo é evitar o esgotamento do storage no `/storage3` sem comprometer a performance de futuras execuções.

---

### 📊 Visão Geral dos Diretórios

| Diretório | Função | Retenção Recomendada |
| :--- | :--- | :--- |
| `singularity_cache/` | Cache das imagens `.sif` (containers). | **Permanente** (Não apagar) |
| `.nextflow_home/` | Plugins, assets e histórico do Nextflow. | **Permanente** (Não apagar) |
| `isugifNF/` | Código-fonte do pipeline clonado. | **Permanente** |
| `nf_work_isugif/` | Arquivos intermediários de cada processo. | **Temporária** (Apagar após validação) |
| `tmp_isugif/` | Swap de I/O, cache de extração e Java TMP. | **Lixo** (Apagar após cada Job) |

---

### 🧹 Protocolos de Limpeza

### 1. Limpeza de Temporários (Imediata)
O diretório definido em `$NXF_TEMP` e `$TMPDIR` acumula resíduos de escrita do Picard/GATK e extrações do Apptainer/Singularity que podem não ser limpos automaticamente se o job sofrer interrupções.

**Comando:**
```bash
rm -rf /storage3/jpitta/oncogensus/tmp_isugif/*
```
### 2. Gerenciamento do Work Directory (-work-dir)
O diretório nf_work_isugif/ é o que mais consome espaço, pois armazena BAMS, FastQs e VCFs intermediários de cada etapa do GATK.

Durante o projeto: Mantenha os arquivos para permitir o uso da flag -resume.

Pós-conclusão: Após mover os resultados finais para uma pasta de entrega/backup, limpe o diretório de trabalho.

Comando manual:

TESTE


