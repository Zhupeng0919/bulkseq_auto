# ============================================================
# Bulkseq 自动化分析流程 — Snakemake 主流程 v2
# ============================================================
# 运行:
#   snakemake --configfile /path/to/project/config.yaml --cores N
# 或:
#   bash run_pipeline.sh /path/to/project/
# ============================================================

import os
import csv

# ---- 从 configfile 加载配置 ----
PROJ  = config["project_dir"]
if not PROJ:
    raise ValueError("config.yaml 中 project_dir 未设置，请填写项目根目录的绝对路径")
PID   = config["project_id"]
MODE  = config.get("mode", "fastq")
SPECIES = config.get("species", "human")

_PIPELINE_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))
_SCRIPTS_DIR  = os.path.join(_PIPELINE_DIR, "scripts")
_CONDA_BIN    = os.path.join(os.environ.get("CONDA_PREFIX", os.path.expanduser("~/miniconda3/envs/bioinfo_env")), "bin")
_RSCRIPT      = os.path.join(_CONDA_BIN, "Rscript")
_LOGS_DIR     = os.path.join(PROJ, "logs")
os.makedirs(_LOGS_DIR, exist_ok = True)

# 参考路径
IDX = ""
GTF = ""
REF_ROOT = config.get("reference", {}).get("root", "")
_SPECIES_TO_GENOME = {"human": "GRCh38.p14", "mouse": "GRCm39", "rat": "GRCr8"}
_GENOME_NAME = _SPECIES_TO_GENOME.get(SPECIES, "")

if REF_ROOT and _GENOME_NAME:
    IDX = config.get("reference", {}).get("index") or os.path.join(REF_ROOT, _GENOME_NAME, f"{_GENOME_NAME}_index")
    GTF = config.get("reference", {}).get("gtf") or os.path.join(REF_ROOT, _GENOME_NAME, f"{_GENOME_NAME}.gtf")
else:
    IDX = config.get("reference", {}).get("index", "")
    GTF = config.get("reference", {}).get("gtf", "")

# 参数
LFC  = config.get("deg", {}).get("logFC_cutoff", 0.5)
PVAL = config.get("deg", {}).get("p_cutoff", 0.05)

DO_GO_KEGG  = config.get("enrichment", {}).get("go_kegg", True)
DO_GSEA     = config.get("enrichment", {}).get("gsea", False)
GSEA_GENES  = config.get("enrichment", {}).get("gsea_target_genes", [])

QC_CHECKPOINT     = config.get("checkpoint", {}).get("after_qc", True)
ALIGN_CHECKPOINT  = config.get("checkpoint", {}).get("after_alignment", True)
SKIP_QC           = config.get("upstream", {}).get("skip_qc", False)

# v2 新增
EXCLUDE_SAMPLES = config.get("exclude_samples", [])
COMPARISON_MODE = config.get("comparison_mode", "auto")
MANUAL_COMPARISONS = config.get("comparisons", [])
# 批次校正
BATCH_CORRECT_ENABLED = config.get("batch", {}).get("correct", False)


# ============================================================
# 辅助函数
# ============================================================

def find_input(in_candidates):
    """在多个候选路径中查找第一个存在的文件"""
    for p in in_candidates:
        expanded = p.format(PROJ = PROJ, PID = PID)
        if os.path.exists(expanded):
            return expanded
    raise FileNotFoundError(f"未找到: {in_candidates}")


def find_counts(wildcards):
    if MODE == "fastq":
        return os.path.join(PROJ, "04_Counts", "counts.txt")
    return find_input(["{PROJ}/04_Counts/counts.txt", "{PROJ}/counts.txt"])


def find_sample_sheet(wildcards):
    return find_input(["{PROJ}/sample_sheet.csv", "{PROJ}/config/sample_sheet.csv"])


def _read_sample_sheet():
    """在 DAG 构建时尝试读取 sample_sheet 以确定分组"""
    candidates = [
        os.path.join(PROJ, "sample_sheet.csv"),
        os.path.join(PROJ, "config", "sample_sheet.csv")
    ]
    for p in candidates:
        if os.path.exists(p):
            with open(p, 'r') as f:
                reader = csv.DictReader(f)
                return list(reader)
    return []


def _read_groups():
    """从 sample_sheet 中提取有序分组列表"""
    rows = _read_sample_sheet()
    groups = []
    for r in rows:
        g = r.get('group', '').strip()
        if g and g not in groups:
            groups.append(g)
    return groups


def _has_batch_column():
    """检测 sample_sheet 是否包含 batch 列"""
    rows = _read_sample_sheet()
    if not rows:
        return False
    return 'batch' in rows[0]


_HAS_BATCH_COL = _has_batch_column()
if BATCH_CORRECT_ENABLED and not _HAS_BATCH_COL:
    raise ValueError(
        "config.yaml 中 batch.correct 设置为 true，"
        "但 sample_sheet.csv 中缺少 'batch' 列。"
        "请在 sample_sheet.csv 中添加 batch 列，或将 batch.correct 设为 false。"
    )
BATCH_CORRECT = BATCH_CORRECT_ENABLED and _HAS_BATCH_COL


def get_comparisons():
    """根据分组和 comparison_mode 返回比较列表。手动 comparions 优先。"""
    # 手动指定的比较列表优先
    if MANUAL_COMPARISONS:
        return list(MANUAL_COMPARISONS)
    groups = _read_groups()
    if len(groups) < 2:
        return []
    mode = COMPARISON_MODE
    if mode == "auto":
        mode = "pairwise" if len(groups) == 2 else "control_vs_rest"

    comps = []
    if mode == "pairwise":
        comps = [f"{groups[1]}_vs_{groups[0]}"]
    elif mode == "control_vs_rest":
        control = groups[0]
        for g in groups[1:]:
            comps.append(f"{g}_vs_{control}")
        if len(groups) > 2:
            comps.append(f"All_vs_{control}")
    elif mode == "all_vs_all":
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                comps.append(f"{groups[j]}_vs_{groups[i]}")
    return comps


# 在 DAG 构建时确定比较列表
COMPARISONS = get_comparisons()


# ============================================================
# 总目标
# ============================================================
rule all:
    input:
        os.path.join(PROJ, "Report.pdf")


# ============================================================
# 模块一：定量数据转换 (04_Counts)
# ============================================================
rule counts_convert:
    input:
        counts = find_counts
    output:
        counts_xls = os.path.join(PROJ, "04_Counts", f"{PID}_counts.xls"),
        fpkm_xls   = os.path.join(PROJ, "04_Counts", f"{PID}_fpkm.xls"),
        tpm_xls    = os.path.join(PROJ, "04_Counts", f"{PID}_tpm.xls")
    params:
        outdir = os.path.join(PROJ, "04_Counts"),
        sdir   = _SCRIPTS_DIR
    shell:
        """
        mkdir -p {params.outdir}
        {_RSCRIPT} {params.sdir}/counts_to_tpm_fpkm.R \
          --input {input.counts} \
          --outdir {params.outdir} \
          --project {PID}
        """


# ============================================================
# 模块 1.5：样本剔除 (04_Counts)
# ============================================================
rule filter_samples:
    """根据 exclude_samples 配置剔除离群样本"""
    input:
        counts = os.path.join(PROJ, "04_Counts", f"{PID}_counts.xls"),
        fpkm   = os.path.join(PROJ, "04_Counts", f"{PID}_fpkm.xls"),
        tpm    = os.path.join(PROJ, "04_Counts", f"{PID}_tpm.xls"),
        sample = find_sample_sheet
    output:
        counts_filt = os.path.join(PROJ, "04_Counts", "filtered", f"{PID}_counts_filtered.xls"),
        fpkm_filt   = os.path.join(PROJ, "04_Counts", "filtered", f"{PID}_fpkm_filtered.xls"),
        tpm_filt    = os.path.join(PROJ, "04_Counts", "filtered", f"{PID}_tpm_filtered.xls"),
        sample_filt = os.path.join(PROJ, "04_Counts", "filtered", "sample_sheet_filtered.csv")
    params:
        outdir  = os.path.join(PROJ, "04_Counts"),
        sdir    = _SCRIPTS_DIR,
        exclude = ",".join(EXCLUDE_SAMPLES) if EXCLUDE_SAMPLES else ""
    shell:
        """
        mkdir -p {params.outdir}/filtered
        {_RSCRIPT} {params.sdir}/filter_samples.R \
          --counts  {input.counts} \
          --sample  {input.sample} \
          --outdir  {params.outdir} \
          --project {PID} \
          --exclude "{params.exclude}"
        """


# ============================================================
# 模块 1.6：批次效应校正 (04_Counts)
# ============================================================
if BATCH_CORRECT:
    rule batch_correct:
        input:
            tpm    = os.path.join(PROJ, "04_Counts", "filtered", f"{PID}_tpm_filtered.xls"),
            sample = os.path.join(PROJ, "04_Counts", "filtered", "sample_sheet_filtered.csv")
        output:
            tpm_bc = os.path.join(PROJ, "04_Counts", f"{PID}_tpm_batch_corrected.xls"),
            done   = touch(os.path.join(PROJ, "04_Counts", ".batch_correct_done"))
        params:
            outdir = os.path.join(PROJ, "04_Counts"),
            sdir   = _SCRIPTS_DIR
        shell:
            """
            {_RSCRIPT} {params.sdir}/batch_correct.R \
              --tpm    {input.tpm} \
              --sample {input.sample} \
              --outdir {params.outdir} \
              --project {PID}
            """


def _deg_inputs(wildcards):
    """返回 deg 规则的输入，当启用批次校正时添加依赖"""
    d = {
        'counts': os.path.join(PROJ, "04_Counts", "filtered", f"{PID}_counts_filtered.xls"),
        'sample': os.path.join(PROJ, "04_Counts", "filtered", "sample_sheet_filtered.csv"),
    }
    if BATCH_CORRECT:
        d['batch_done'] = os.path.join(PROJ, "04_Counts", ".batch_correct_done")
    return d


# ============================================================
# 模块二：差异表达分析 (05_DEG)
# ============================================================
rule deg:
    input:
        unpack(_deg_inputs)
    output:
        done    = os.path.join(PROJ, "05_DEG", ".deg_done"),
        alldiff = os.path.join(PROJ, "05_DEG", f"{PID}_All_Diff.csv")
    params:
        outdir    = os.path.join(PROJ, "05_DEG"),
        sdir      = _SCRIPTS_DIR,
        comp_mode = COMPARISON_MODE,
        comps     = ",".join(COMPARISONS)
    shell:
        """
        mkdir -p {params.outdir}
        {_RSCRIPT} {params.sdir}/deg_deseq2.R \
          --counts   {input.counts} \
          --sample   {input.sample} \
          --outdir   {params.outdir} \
          --project  {PID} \
          --lfc {LFC} --pval {PVAL} \
          --comparison_mode {params.comp_mode} \
          --comparisons "{params.comps}"
        touch {output}
        """


# ============================================================
# 模块三：GO / KEGG 富集 (06_GO_KEGG)
# ============================================================
if DO_GO_KEGG and len(COMPARISONS) > 0:
    rule go_kegg:
        input:
            deg_done = os.path.join(PROJ, "05_DEG", ".deg_done")
        output:
            touch(os.path.join(PROJ, "06_GO_KEGG", ".gokegg_done"))
        params:
            outdir = os.path.join(PROJ, "06_GO_KEGG"),
            comps  = COMPARISONS
        run:
            os.makedirs(params.outdir, exist_ok = True)
            for comp in params.comps:
                diff_csv = os.path.join(PROJ, "05_DEG", f"{PID}_{comp}_Diff.csv")
                shell("""
                    {{_RSCRIPT}} {0}/go_kegg.R \
                      --deg        {1} \
                      --outdir     {2} \
                      --project    {3} \
                      --species    {4} \
                      --comparison {5} \
                      --pval       {6}
                """.format(_SCRIPTS_DIR, diff_csv, params.outdir, PID, SPECIES, comp, PVAL))
else:
    rule go_kegg:
        input:
            deg_done = os.path.join(PROJ, "05_DEG", ".deg_done")
        output:
            touch(os.path.join(PROJ, "06_GO_KEGG", ".gokegg_done"))
        run:
            os.makedirs(os.path.join(PROJ, "06_GO_KEGG"), exist_ok = True)


# ============================================================
# 模块四：GSEA 分析 (07_GSEA)
# ============================================================
if DO_GSEA and len(COMPARISONS) > 0:
    rule gsea:
        input:
            deg_done = os.path.join(PROJ, "05_DEG", ".deg_done")
        output:
            touch(os.path.join(PROJ, "07_GSEA", ".gsea_done"))
        params:
            outdir = os.path.join(PROJ, "07_GSEA"),
            comps  = COMPARISONS
        run:
            os.makedirs(params.outdir, exist_ok = True)
            for comp in params.comps:
                diff_csv = os.path.join(PROJ, "05_DEG", f"{PID}_{comp}_Diff.csv")
                shell("""
                    {{_RSCRIPT}} {0}/gsea.R \
                      --deg        {1} \
                      --outdir     {2} \
                      --project    {3} \
                      --species    {4} \
                      --comparison {5}
                """.format(_SCRIPTS_DIR, diff_csv, params.outdir, PID, SPECIES, comp))
else:
    rule gsea:
        input:
            deg_done = os.path.join(PROJ, "05_DEG", ".deg_done")
        output:
            touch(os.path.join(PROJ, "07_GSEA", ".gsea_done"))
        run:
            os.makedirs(os.path.join(PROJ, "07_GSEA"), exist_ok = True)
            shell("touch {0}/SKIPPED".format(os.path.join(PROJ, "07_GSEA")))


# ============================================================
# 模块五：PDF 报告生成
# ============================================================
def _report_inputs(wildcards):
    """返回 report 规则的输入，优先使用批次校正后的 TPM"""
    d = {
        'deg_done':  os.path.join(PROJ, "05_DEG", ".deg_done"),
        'gokegg':    os.path.join(PROJ, "06_GO_KEGG", ".gokegg_done"),
        'gsea_done': os.path.join(PROJ, "07_GSEA", ".gsea_done"),
    }
    if BATCH_CORRECT:
        d['tpm'] = os.path.join(PROJ, "04_Counts", f"{PID}_tpm_batch_corrected.xls")
    else:
        d['tpm'] = os.path.join(PROJ, "04_Counts", "filtered", f"{PID}_tpm_filtered.xls")
    return d


rule report:
    input:
        unpack(_report_inputs)
    output:
        pdf = os.path.join(PROJ, "Report.pdf")
    params:
        rmd         = os.path.join(_SCRIPTS_DIR, "generate_report.Rmd"),
        comparisons = ",".join(COMPARISONS),
        comp_mode   = COMPARISON_MODE,
        exclude     = ",".join(EXCLUDE_SAMPLES) if EXCLUDE_SAMPLES else ""
    shell:
        """
        export PATH="{_CONDA_BIN}:$PATH"
        {_RSCRIPT} -e 'rmarkdown::render("{params.rmd}",
          params = list(
            project_dir     = "{PROJ}",
            project_id      = "{PID}",
            mode            = "{MODE}",
            species         = "{SPECIES}",
            comparisons     = "{params.comparisons}",
            comparison_mode = "{params.comp_mode}",
            exclude_samples = "{params.exclude}"
          ),
          output_file = "Report.pdf",
          output_dir  = "{PROJ}",
          quiet = TRUE
        )'
        rm -f {PROJ}/04_Counts/counts.txt
        rm -f {PROJ}/04_Counts/counts.txt.summary
        rm -rf {_LOGS_DIR}
        """


# ============================================================
# MODE 1 RULES (FASTQ → Counts)  — 上游流程
# ============================================================
if MODE == "fastq":

    RAW   = os.path.join(PROJ, "RawData")
    CLEAN = os.path.join(PROJ, "01_CleanData")
    QC    = os.path.join(PROJ, "02_QC_Reports")
    ALIGN = os.path.join(PROJ, "03_Alignment")
    COUNT = os.path.join(PROJ, "04_Counts")

    if not IDX:
        raise ValueError(
            "mode=fastq 需要 HISAT2 索引，但 config.yaml 中 reference.index 为空。\n"
            "请先构建参考基因组: bash scripts/build_reference.sh <species>")
    if not GTF:
        raise ValueError(
            "mode=fastq 需要基因注释 GTF，但 config.yaml 中 reference.gtf 为空。\n"
            "请先构建参考基因组: bash scripts/build_reference.sh <species>")

    def get_samples():
        """
        扫描 RawData/ 下的 PE FASTQ 文件，返回样本列表。
        支持两种布局：
          - 子目录布局：RawData/<样本名>/*_1.fq.gz （子目录名即为样本名）
          - 扁平布局：  RawData/<样本名>_1.fq.gz （文件名前缀为样本名）
        """
        import glob
        samples = set()
        raw_dir = os.path.join(PROJ, "RawData")
        if not os.path.isdir(raw_dir):
            return []

        # 优先扫描子目录布局
        for entry in sorted(os.listdir(raw_dir)):
            entry_path = os.path.join(raw_dir, entry)
            if os.path.isdir(entry_path):
                r1_files = glob.glob(os.path.join(entry_path, "*_1.fq.gz"))
                if r1_files:
                    samples.add(entry)

        # 回退：扁平布局
        if not samples:
            for fq in glob.glob(os.path.join(raw_dir, "*_1.fq.gz")):
                base = os.path.basename(fq)
                samples.add(base.replace("_1.fq.gz", ""))

        return sorted(samples)

    def get_fastq_r1(wildcards):
        """解析样本的 R1 FASTQ 路径。先查子目录，再查扁平路径。"""
        import glob
        sample = wildcards.sample
        subdir = os.path.join(PROJ, "RawData", sample, "*_1.fq.gz")
        flat   = os.path.join(PROJ, "RawData", f"{sample}_1.fq.gz")
        files  = glob.glob(subdir)
        if files:
            return sorted(files)[0]
        if os.path.exists(flat):
            return flat
        raise FileNotFoundError(f"未找到样本 '{sample}' 的 R1 FASTQ 文件")

    def get_fastq_r2(wildcards):
        """由 R1 路径推导 R2 路径。"""
        r1 = get_fastq_r1(wildcards)
        r2 = r1.replace("_1.fq.gz", "_2.fq.gz")
        if os.path.exists(r2):
            return r2
        raise FileNotFoundError(f"未找到配对 R2 文件: {r2}")

    rule fastp:
        input:
            r1 = get_fastq_r1,
            r2 = get_fastq_r2
        output:
            r1c = os.path.join(CLEAN, "{sample}_1.clean.fq.gz"),
            r2c = os.path.join(CLEAN, "{sample}_2.clean.fq.gz"),
            html = os.path.join(QC, "{sample}_fastp.html")
        params:
            trim = config.get("upstream", {}).get("trim_front", 15)
        shell:
            """
            mkdir -p {CLEAN} {QC}
            fastp -i {input.r1} -I {input.r2} \
                  -o {output.r1c} -O {output.r2c} \
                  --trim_front1 {params.trim} --trim_front2 {params.trim} \
                  -w 4 \
                  -h {output.html} -j {QC}/{wildcards.sample}.json
            """

    rule fastqc:
        input:
            r1 = os.path.join(CLEAN, "{sample}_1.clean.fq.gz"),
            r2 = os.path.join(CLEAN, "{sample}_2.clean.fq.gz")
        output:
            r1qc = os.path.join(QC, "{sample}_1.clean_fastqc.html"),
            r2qc = os.path.join(QC, "{sample}_2.clean_fastqc.html")
        shell:
            """
            fastqc {input.r1} {input.r2} -t 2 -o {QC}
            """

    rule multiqc:
        input:
            expand(os.path.join(QC, "{s}_fastp.html"), s = get_samples()),
            expand(os.path.join(QC, "{s}_1.clean_fastqc.html"), s = get_samples()),
            expand(os.path.join(QC, "{s}_2.clean_fastqc.html"), s = get_samples())
        output:
            report = os.path.join(QC, "multiqc_final_report.html")
        shell:
            """
            multiqc {QC} -o {QC} -n multiqc_final_report --force
            """

    rule qc_checkpoint:
        input:
            os.path.join(QC, "multiqc_final_report.html")
        output:
            touch(os.path.join(QC, ".qc_done"))
        run:
            if QC_CHECKPOINT:
                print("\n" + "=" * 70)
                print("  CHECKPOINT 1: QC 质量审查")
                print("=" * 70)
                print(f"  MultiQC 报告: {QC}/multiqc_final_report.html")
                print()
                print("  请检查 MultiQC 报告中的以下指标：")
                print("    - Per base sequence quality (每个碱基的质量分数)")
                print("    - Per sequence quality scores (每条序列的质量)")
                print("    - Adapter content (接头污染比例)")
                print("    - GC content (GC 含量分布)")
                print("    - Sequence duplication levels (序列重复水平)")
                print()
                print("  如果 QC 质量较差（大量低质量碱基/高接头污染），")
                print("  可以在后续步骤中启用多比对模式以保留更多 reads：")
                print("    在 config.yaml 中设置:")
                print("      upstream:")
                print("        multimapping: true")
                print()
                print("  确认 QC 质量合格后，将 config.yaml 中的:")
                print("      checkpoint:")
                print("        after_qc: false")
                print("  并重新运行 pipeline。")
                print("=" * 70 + "\n")
                import sys
                sys.exit(1)

    rule hisat2:
        input:
            r1 = os.path.join(CLEAN, "{sample}_1.clean.fq.gz"),
            r2 = os.path.join(CLEAN, "{sample}_2.clean.fq.gz"),
            qc_ok = os.path.join(QC, ".qc_done")
        output:
            bam = os.path.join(ALIGN, "{sample}.sorted.bam"),
            bai = os.path.join(ALIGN, "{sample}.sorted.bam.bai"),
            log = os.path.join(_LOGS_DIR, "{sample}_hisat2.log")
        params:
            index = IDX
        shell:
            """
            mkdir -p {ALIGN} {_LOGS_DIR}
            hisat2 -p 10 -x {params.index} \
                   -1 {input.r1} -2 {input.r2} \
                   2> {output.log} | \
            samtools view -@ 2 -bS - | \
            samtools sort -@ 4 -m 4G -o {output.bam} -
            samtools index {output.bam}
            """

    rule mapping_summary:
        input:
            logs = expand(os.path.join(_LOGS_DIR, "{s}_hisat2.log"), s = get_samples())
        output:
            summary = os.path.join(PROJ, "total_mapping_summary.txt")
        shell:
            """
            echo -e "SampleID\\tHISAT2_Mapping_Rate" > {output.summary}
            for log in {_LOGS_DIR}/*_hisat2.log; do
              s_name=$(basename "$log" _hisat2.log)
              h_rate=$(grep "overall alignment rate" "$log" | awk '{{print $1}}')
              echo -e "${{s_name}}\\t${{h_rate}}" >> {output.summary}
            done
            """

    rule alignment_checkpoint:
        input:
            os.path.join(PROJ, "total_mapping_summary.txt")
        output:
            touch(os.path.join(ALIGN, ".align_done"))
        run:
            if ALIGN_CHECKPOINT:
                summary_file = os.path.join(PROJ, "total_mapping_summary.txt")
                print("\n" + "=" * 70)
                print("  CHECKPOINT 2: 比对质量审查")
                print("=" * 70)
                print(f"  比对统计文件: {summary_file}")
                print(f"  HISAT2 比对日志: {_LOGS_DIR}/")
                print()
                print("  请检查以下内容：")
                print("    1. 查看 total_mapping_summary.txt 中各样本的比对率")
                print("    2. 正常 RNA-seq 比对率应在 70%-95% 范围")
                print("    3. 如果某个样本比对率异常低 (< 50%)，检查 hisat2 日志")
                print(f"       {_LOGS_DIR}/<样本名>_hisat2.log")
                print()
                if os.path.exists(summary_file):
                    print("  当前比对率摘要：")
                    print("  " + "-" * 50)
                    with open(summary_file) as f:
                        for line in f:
                            print(f"  {line.rstrip()}")
                    print("  " + "-" * 50)
                    print()
                print("  确认比对质量合格后，将 config.yaml 中的:")
                print("      checkpoint:")
                print("        after_alignment: false")
                print("  并重新运行 pipeline。")
                print("=" * 70 + "\n")
                import sys
                sys.exit(1)

    rule featurecounts:
        input:
            bams = expand(os.path.join(ALIGN, "{s}.sorted.bam"), s = get_samples()),
            summary_rule = os.path.join(PROJ, "total_mapping_summary.txt"),
            align_ok = os.path.join(ALIGN, ".align_done")
        output:
            counts = os.path.join(COUNT, "counts.txt")
        params:
            gtf  = GTF,
            log  = os.path.join(_LOGS_DIR, "fc.log"),
            multi = "" if not config.get("upstream", {}).get("multimapping") else "-M --fraction"
        shell:
            """
            mkdir -p {COUNT} {_LOGS_DIR}
            featureCounts -T 16 -p -B -C -s 0 {params.multi} \
                          -a {params.gtf} \
                          -o {output.counts} {ALIGN}/*.sorted.bam \
                          2> {params.log}
            """
