# WGS de caprinos, bovinos e ovinos — guia de execução

Código em inglês; instruções em português. Pipeline para reads curtas Illumina paired-end de animais modernos. Não é um pipeline de RNA-seq, reads longas ou DNA antigo. Vai dos FASTQs ao BAM; não faz chamada de variantes.

## 1. O que mudar

Usar uma diretoria independente por espécie/projeto. Não misturar referências numa execução.

| Ficheiro | Campo | Alteração |
|---|---|---|
| settings.env | SPECIES | goat, cattle ou sheep; etiqueta informativa, não seleciona referência automaticamente |
| settings.env | REFERENCE | caminho absoluto para o FASTA da montagem escolhida para a espécie |
| settings.env | SLURM_ACCOUNT | conta SLURM do colega; não presumir que é a conta de Fábio |
| settings.env | SLURM_PARTITION | medium, ou outra partição autorizada com recursos adequados |
| settings.env | NEXTFLOW_ACTIVATE | script de ativação Nextflow acessível pelo colega |
| settings.env | ADAPTERS_BASE | FASTA de adaptadores correspondente à preparação das bibliotecas |
| samples.tsv | todas as linhas | animais, runs, bibliotecas e caminhos dos FASTQs |
| nextflow.config | caminhos de ferramentas | manter se o colega tiver acesso; corrigir se instalações pessoais não forem legíveis/executáveis |

**Não existe parâmetro de número de cromossomas neste pipeline.** O BWA usa as sequências do FASTA, e o SAMtools usa o cabeçalho BAM. Não é preciso adaptar `--chr-set`: essa opção seria relevante noutras ferramentas, como PLINK, em análises posteriores. Todas as sequências da referência entram na profundidade média, incluindo scaffolds e cromossomas sexuais. Não é uma média apenas dos autossomas.

Escolher e registar a montagem e a versão exatas; preferir uma referência consistente com as análises futuras. O pipeline não escolhe nem descarrega automaticamente uma montagem. FASTA não comprimido, com índices BWA correspondentes e cabeçalhos sem identificadores duplicados.

## 2. Instalação e módulos

Extrair no HPC, numa pasta própria:

```bash
tar -xzf livestock_wgs_nextflow.tar.gz
cd livestock_wgs_nextflow
cp settings.env.example settings.env
nano settings.env
```

Mesmo HPC: preservados Java/21.0.8, GCC/14.3.0, SAMtools/1.22.1, BCFtools/1.22 e ativação do ambiente Conda BWA através de site/langs e Anaconda3/2025.12-1. FastQC é executável local; Trimmomatic e Picard são JARs. Ver caminhos em nextflow.config. As versões efetivas são registadas nos outputs das ferramentas quando disponíveis.

Os caminhos das ferramentas são exemplos e devem ser configurados para o cluster de destino. Estar no mesmo HPC NÃO garante permissões nessas pastas: confirmar acesso nos nós de cálculo. Se necessário, usar instalações próprias com as mesmas versões e atualizar os caminhos. Não é preciso instalar módulos diferentes por espécie.

## 3. Preparar os índices da referência (uma vez)

Se já existirem `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` junto do FASTA, não repetir a indexação. Caso contrário, submeter um job com memória e tempo suficientes:

```bash
source settings.env
sbatch --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" \
  --cpus-per-task=2 --mem=32G --time=12:00:00 \
  --output=reference_index_%j.log index_reference.sh "$REFERENCE"
```

O FASTA deve estar numa diretoria onde o utilizador tenha permissão para escrever os índices. O script também cria `.fai` com samtools faidx. Não executar duas indexações simultâneas sobre a mesma referência.

## 4. Preparar samples.tsv

```bash
cp samples.example.tsv samples.tsv
nano samples.tsv
```

É TSV: separar colunas com TAB, não espaços. Uma linha por run, não uma linha por mate.

- `animal_id`: identificador biológico único; repetir nas linhas dos runs do mesmo animal.
- `run`: identificador único para cada run; usado nos nomes dos ficheiros e read group ID.
- `library_id`: mesma etiqueta APENAS para runs da mesma biblioteca física. Usar etiquetas diferentes para bibliotecas distintas. Se desconhecida, deixar vazio: usa o run como biblioteca e não assume duplicação entre runs. Preferir nomes únicos por animal, por exemplo ANIMAL01_LIB1.
- `header_format`: normalmente `standard`. `hash_mate` apenas para o formato especial `@...#1:Y:0:...` / `#2:...`; cria derivados normalizados e preserva os originais.
- `R1`, `R2`: caminhos absolutos dos FASTQs correspondentes, previamente verificados e sincronizados.

Não repetir o mesmo FASTQ em dois runs. Usar letras, números, `_`, `-` ou `.` nos identificadores, sem espaços. O manifest de exemplo tem dois animais, um com dois runs; substituir todas as linhas fictícias.

## 5. Executar

Na pasta do pipeline:

```bash
nohup bash run.sh "$PWD/samples.tsv" --stage full -resume \
  > run_full.log 2>&1 < /dev/null &
```

O controlador Nextflow submete as tarefas ao SLURM; não submeter cada processo à mão. Seguir as regras locais sobre onde manter o controlador. Não iniciar dois controladores na mesma diretoria.

```bash
tail -n 40 run_full.log
squeue -u "$USER"
```

Outras etapas possíveis:

```bash
bash run.sh "$PWD/samples.tsv" --stage raw_qc -resume
bash run.sh "$PWD/samples.tsv" --stage trim_qc -resume
```

`raw_qc` termina após FastQC raw; `trim_qc` termina após trimming e QC; `full` continua até BAM e estatísticas. Ao passar para full na mesma pasta, `-resume` reutiliza tarefas concluídas cuja cache continue válida. Não mover inputs nem apagar work/. São até 15 tarefas simultâneas, não 15 animais completos. Alterar `executor.queueSize` em nextflow.config se necessário; verificar limites de CPU/memória da conta.

## 6. QC e revisão manual

Fluxo: FASTQC_RAW → decisão → normalização opcional → TRIM → FASTQC_TRIM → decisão → MAP → merge de todos os runs do animal → Picard → estatísticas.

- RAW: apenas FAIL em Adapter Content retém o run, a menos que seja autorizado após revisão.
- Após trimming: FAIL em Adapter Content OU Per base sequence quality num dos mates retém o run inteiro.
- WARN e todos os outros módulos deixam passar. Esta é uma política do projeto, não uma garantia de ausência de outros problemas.
- Não há bypass do controlo pós-trimming. Falta de módulos, contagens incompatíveis ou ausência de pares também impedem mapping.
- A retenção por QC não interrompe os restantes runs. Erros de execução/software continuam a interromper o workflow e exigem correção.
- O merge de um animal só avança com TODOS os runs previstos no manifest; os outros animais continuam independentemente.

Ficheiros editáveis, após terminar a execução:

1. `additional_adapters.fa`: inicialmente vazio. Acrescentar adaptadores confirmados em FASTA; não copiar automaticamente qualquer sequência sobre-representada (pode ser biológica).
2. `approved_raw_adapter_runs.txt`: um run por linha, autorizado para trimming após inspeção do FAIL raw. Também permite usar adaptadores adicionais num run com raw aprovado.
3. `trim_overrides.tsv`: uma linha por run com opções próprias. Colunas `run` e `options`, separadas por TAB.

Trimming padrão para TODAS as espécies:

```text
ILLUMINACLIP:<adapters>:2:30:10:1:true AVGQUAL:20 MINLEN:50
```

Não há tratamentos especiais herdados de amostras de suínos. Para mudar um run, colocar por exemplo em `options`:

```text
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```

As opções substituem a parte de qualidade/comprimento, mantendo ILLUMINACLIP. Para alterar o padrão de todos os runs, editar `params.quality_options` em nextflow.config. A escolha deve depender das bibliotecas e do QC, não apenas da espécie.

A lista adicional é partilhada apenas pelos runs aprovados na lista; os outros usam só os adaptadores padrão. Editar a lista adicional e retomar pode repetir trimming dos runs que a usam, mas não dos restantes. Só avançam ao mapping se passarem o QC pós-trimming. As listas e o FASTA usado ficam registados em provenance/.

## 7. BAMs e duplicados

MAP usa BWA-MEM, apenas os pares sobreviventes; unpaired são guardados mas não mapeados. Read groups: SM=animal_id, ID=run, LB=library_id (ou run quando não indicado). Ordenação/indexação: SAMtools. Runs da mesma amostra são fundidos antes de Picard.

Picard: `REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true`. Duplicados de biblioteca/PCR são marcados e permanecem; os identificados como duplicados de sequenciação/ópticos são removidos. `.dedup.bam` NÃO significa que todos os duplicados PCR foram removidos. Java usa heap de 24 GB, com job de 64 GB, conforme o ajuste feito no HPC.

## 8. Resultados e relatórios

- `results/01_fastqc_raw/`: HTML/ZIP raw por run.
- `results/02_trim/`: paired/unpaired, log de trimming e contagens.
- `results/03_fastqc_trimmed/`: HTML/ZIP após trimming.
- `results/04_mapping/`: BAM ordenado por run.
- `results/04_merged/`: BAM por animal antes do Picard.
- `results/05_duplicates/`: BAM final, índice, métricas e log Picard.
- `results/stats/`: detalhes por run, merged e final.
- `results/tables/run_statistics.tsv`: contagens e profundidade por run.
- `results/tables/animal_statistics.tsv`: runs esperados/concluídos/pendentes, mapped reads e profundidade antes/depois do Picard.
- `results/tables/qc_review.tsv`: flags e estado dos runs.
- `provenance/<execução>/report.html`, `timeline.html`, `trace.tsv`: execução/recursos do Nextflow, não as profundidades biológicas.

Contagens de reads são reads individuais, não pares. `reads_after_trim` inclui unpaired; `reads_sent_to_mapping` inclui só paired. Mapped reads exclui alinhamentos secundários/suplementares e não impõe limiar MAPQ.

Profundidade: `samtools depth -aa -s -q 20 -Q 20 -G 3844`, inclui zeros em toda a referência, qualidade de base e mapeamento >=20, sobreposição do par contada uma vez; exclui flags de duplicado, QC fail, secundários/suplementares e unmapped. Antes do Picard as flags de duplicados ainda não foram atribuídas. `final_mean_depth` é a média após Picard. A média conjunta é recalculada no BAM, não é a média aritmética das médias por run.

Tabelas agregadas aparecem ao concluir a execução elegível. Estado FINAL indica conclusão das estatísticas do animal. Uma execução pode concluir e deixar runs retidos. Não inferir aprovação atual só pela existência de um BAM de execução anterior.

Outputs grandes são links para work/: NÃO apagar work/ nem usar nextflow clean até copiar e verificar os resultados finais. Não há eliminação automática dos dados de entrada.

## 9. Retomar e partilhar

Depois de corrigir uma falha ou os parâmetros de um run, executar novamente o comando com `-resume`. Só tarefas com cache válida são reaproveitadas; alterações de comandos/inputs podem invalidar etapas e descendentes.

Para GitHub, incluir código, configurações de exemplo, guia e testes; não incluir FASTQs/BAMs, work/, results/ ou provenance/. Rever settings.env antes de publicar: os caminhos e conta são específicos do utilizador. O .gitignore exclui outputs e caches. Não foi publicado nenhum repositório por este pacote.

## 10. Validação e fontes

```bash
python3 -m unittest discover -s tests -v
```

Testes locais verificam parsing QC, WARN/FAIL, contagens e cálculo de profundidade. O fluxo de mapeamento que serviu de base já foi executado no HPC com suínos. Esta versão genérica ainda não foi executada de ponta a ponta com dados de caprinos/bovinos/ovinos. Começar com um animal pequeno antes de lançar o conjunto completo.

- BWA: https://github.com/lh3/bwa
- SAMtools: https://www.htslib.org/doc/
- Picard: https://broadinstitute.github.io/picard/command-line-overview.html
- Nextflow: https://www.nextflow.io/docs/latest/
