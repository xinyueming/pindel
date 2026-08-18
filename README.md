# Pindel

Pindel 是一款基于 split-read（分裂读）算法的结构变异（Structural Variant, SV）检测工具，专门用于检测：

- **缺失 (Deletions, DEL)**
- **短插入 (Short Insertions, SI)**
- **长插入 (Long Insertions, LI)**
- **倒位 (Inversions, INV)**
- **串联重复 (Tandem Duplications, TD)**
- **移动元件插入 (Mobile Element Insertions, MEI)**

## 安装

### pip 安装（推荐）

```bash
pip install pindel-tool
```

安装后可直接使用 `pindel-tool` 命令：

```bash
pindel-tool --help
pindel-tool pindel --help
pindel-tool filter --help
```

### 源码构建（开发者）

```bash
pip install build
python -m build
pip install dist/pindel_tool-*.whl
```

## 目录结构

```
pindel/
├── pyproject.toml          # Python 包配置
├── README.md               # 本文档
├── COPYING.txt             # 许可证
├── .gitignore
└── src/
    └── pindel_tool/        # Python 包
        ├── __init__.py
        ├── cli.py          # CLI 入口
        ├── pindel          # Pindel 二进制
        ├── pindel2vcf      # VCF 转换二进制
        └── annovar_scripts/
            └── table_annovar.pl
```

## 使用方法

### 1. 准备输入文件

#### 配置文件格式

创建一个配置文件（如 `config.txt`），每行包含 BAM 文件路径、插入片段大小和样本名（以制表符分隔）：

```
/path/to/sample1.bam    350    SAMPLE1
/path/to/sample2.bam    350    SAMPLE2
/path/to/tumor.bam      400    TUMOR
/path/to/normal.bam     400    NORMAL
```

**说明：**
- **BAM 文件**：必须预先用 `samtools index` 建立 `.bai` 索引
- **插入片段大小 (Insert Size)**：配对末端测序的平均插入片段长度，可从测序人员处获取，不确定时可用 500 作为默认值
- **样本名**：输出结果中的样本标识

### 2. 基本运行

```bash
# 基础命令
pindel-tool pindel -f reference.fa -i config.txt -o output_prefix -c ALL

# 参数说明：
# -f/--fasta      参考基因组 FASTA 文件
# -i/--config     配置文件路径
# -o/--output     输出文件前缀
# -c/--chromosome 染色体名称或 ALL（全部染色体）
```

### 3. 常用参数

#### 必需参数

| 参数 | 长参数 | 说明 |
|------|--------|------|
| `-f` | `--fasta` | 参考基因组 FASTA 文件 |
| `-i` | `--config-file` | BAM 配置文件 |
| `-o` | `--output-prefix` | 输出文件前缀 |
| `-c` | `--chromosome` | 染色体名称（如 20）或 ALL |

#### 可选参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-T` | 线程数 | 1 |
| `-j` | 指定分析的染色体区域（BED 格式） | - |
| `-J` | 排除指定染色体区域 | - |
| `-l` | 输出长插入（LI）结果 | 关闭 |
| `-q` | 检测移动元件插入（MEI） | 关闭 |
| `-b` | 使用 BreakDancer 结果辅助检测 | - |
| `-g` | 进行基因分型（genotyping） | 关闭 |
| `-e` | 最小支持读数阈值 | 1 |
| `-W` | 窗口大小（影响内存使用） | 默认 |
| `-a` | 最小删除长度 | 1 |
| `-z` | 最大删除长度 | 默认无限制 |

#### 区域限制参数（节省内存和时间）

```bash
# 分析单个染色体
pindel-tool pindel -f ref.fa -i config.txt -o output -c 20

# 分析染色体特定区域
pindel-tool pindel -f ref.fa -i config.txt -o output -c 20:1000000-2000000

# 排除着丝粒/端粒区域（减少假阳性）
pindel-tool pindel -f ref.fa -i config.txt -o output -c ALL -J centromere.bed
```

### 4. 性能优化

#### 多线程加速

```bash
# 使用 4 线程
pindel-tool pindel -f ref.fa -i config.txt -o output -c ALL -T 4
```

#### 并行处理（集群环境）

```bash
# 按染色体并行提交任务
for chr in {1..22} X Y; do
    pindel-tool pindel -f ref.fa -i config.txt -o chr${chr} -c ${chr} &
done
wait
```

#### 降低内存使用

```bash
# 减小窗口大小
pindel-tool pindel -f ref.fa -i config.txt -o output -c ALL -W 1

# 排除高重复区域
pindel-tool pindel -f ref.fa -i config.txt -o output -c ALL -J exclude.bed
```

## 输出文件

Pindel 输出多个文本文件，每种 SV 类型一个：

| 文件后缀 | SV 类型 | 说明 |
|----------|---------|------|
| `_D` | Deletion | 缺失 |
| `_SI` | Short Insertion | 短插入（<50bp） |
| `_LI` | Long Insertion | 长插入（需 `-l` 参数） |
| `_INV` | Inversion | 倒位 |
| `_TD` | Tandem Duplication | 串联重复 |
| `_BP` | Breakpoint | 断点文件 |
| `_MEI` | Mobile Element Insertion | 移动元件插入（需 `-q` 参数） |

### 输出格式解析

每行包含 SV 事件信息，格式如下：

```
事件序号  SV类型  起始位置  终止位置  长度  样本名  左断点深度  右断点深度  +链支持数  +链独特读数  -链支持数  -链独特读数
```

**字段说明：**
- **左/右断点深度**：该位置的测序深度，若远高于支持读数，可能为假阳性
- **+/-链支持数**：支持该 SV 的正/负链 read 总数
- **独特读数**：去除 PCR 重复后的独立支持读数

## 转换为 VCF 格式

使用 `pindel-tool pindel2vcf` 将 Pindel 输出转换为 VCF：

```bash
pindel-tool pindel2vcf -r reference.fa -R GRCh38 -d 20200101 -P output_D -e 5

# 参数说明：
# -r  参考基因组 FASTA 文件
# -R  参考基因组名称
# -d  参考基因组版本日期
# -P  Pindel 输出文件前缀（如 output 对应 output_D, output_SI 等）
# -e  最小支持读数阈值（过滤低质量结果）
# -v  输出 VCF 文件名（可选，默认使用输入文件名）
# -G  输出 GATK 兼容格式
```

### 示例：完整分析流程

```bash
# 1. 创建配置文件
echo "/data/tumor.bam    400    TUMOR" > config.txt
echo "/data/normal.bam   400    NORMAL" >> config.txt

# 2. 运行 Pindel（排除着丝粒区域）
pindel-tool pindel -f /data/reference.fa -i config.txt -o pindel_result -c ALL -T 4 -J centromere.bed

# 3. 转换为 VCF（过滤支持读数 <5 的结果）
pindel-tool pindel2vcf -r /data/reference.fa -R GRCh38 -d 20200101 -P pindel_result -e 5 -v deletions.vcf
pindel-tool pindel2vcf -r /data/reference.fa -R GRCh38 -d 20200101 -P pindel_result -e 5 -v inversions.vcf
```

## Python CLI 工具

安装后可通过 `pindel-tool` 命令运行完整分析流程：

```bash
pip install pindel-tool
```

### 四个子命令

| 子命令 | 功能 |
|--------|------|
| `pindel` | 调用 pindel 检测结构变异 |
| `pindel2vcf` | 将 pindel 输出转为 VCF |
| `anno` | ANNOVAR 基因注释 |
| `filter` | 按基因/转录本过滤 VCF，输出表格 |

### 完整流程示例

```bash
# 1. 检测结构变异
pindel-tool pindel -f reference.fa -i config.txt -o output -c ALL -T 4

# 2. 转换为 VCF
pindel-tool pindel2vcf -P output -r reference.fa -R hg38 -d 20230801 -v result.vcf

# 3. ANNOVAR 注释
pindel-tool anno result.vcf /path/to/humandb/ --buildver hg38

# 4. 按基因过滤输出
pindel-tool filter result.hg38_multianno.vcf --gene FLT3,KMT2A -o filtered.tsv
```

### 子命令参数

#### pindel

```bash
pindel-tool pindel -f <reference.fa> -i <config.txt> -o <prefix> [options]
```

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-f` | 参考基因组 FASTA | 必需 |
| `-i` | BAM 配置文件 | 必需 |
| `-o` | 输出文件前缀 | 必需 |
| `-c` | 染色体（或 ALL） | ALL |
| `-T` | 线程数 | 1 |
| `-a` | 错配阈值 | 1 |
| `-M` | 最小支持读数 | 1 |

#### pindel2vcf

```bash
pindel-tool pindel2vcf -P <prefix> -r <ref.fa> -R <name> -d <date> [options]
```

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-P` | pindel 输出文件前缀 | 必需 |
| `-r` | 参考基因组 FASTA | 必需 |
| `-R` | 参考基因组名称（如 hg38） | 必需 |
| `-d` | 参考基因组版本日期 | 必需 |
| `-v` | 输出 VCF 文件名 | 自动生成 |
| `-e` | 最小支持读数 | 1 |
| `-he` | 杂合阈值 | 0.2 |
| `-ho` | 纯合阈值 | 0.8 |

#### anno

```bash
pindel-tool anno <input.vcf> [db_path] [options]
```

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `input_vcf` | 输入 VCF 文件 | 必需 |
| `db_path` | ANNOVAR 数据库路径 | `/mnt/nas/zhangrs/3.database/hg38/humandb/` |
| `--buildver` | 基因组版本 | hg38 |
| `--protocol` | 注释协议 | refGeneWithVer |
| `--operation` | 操作类型 | g |
| `--nastring` | 空值标记 | . |
| `--out` | 输出文件前缀 | 基于输入文件名 |

#### filter

```bash
pindel-tool filter <annotated.vcf> [options]
```

| 参数 | 说明 |
|------|------|
| `input_vcf` | 注释后的 VCF（`*_multianno.vcf`） |
| `--gene` | 基因名，逗号分隔（如 `FLT3,KMT2A`） |
| `--transcript` | 转录本 ID，逗号分隔 |
| `--gene-transcript-pair` | 基因:转录本配对（如 `FLT3:NM_004119.3`），可重复指定 |
| `-o` | 输出文件路径（默认 stdout） |

**筛选逻辑**：多种筛选条件之间为 OR 关系。

**输出字段**：CHROM, POS, ID, REF, ALT, Gene, Transcript, SVTYPE, SVLEN, Insertion, CDS, AA, GT, AD, VD, DP, AF, Sample

## 肿瘤-正常配对分析

对于肿瘤样本的体细胞突变检测：

```bash
# 创建配置文件
echo "/data/tumor.bam    350    TUMOR" > tn_config.txt
echo "/data/normal.bam   350    NORMAL" >> tn_config.txt

# 运行 Pindel
pindel-tool pindel -f hg38.fa -i tn_config.txt -o somatic_sv -c ALL -T 4 -g

# 后续需要比较肿瘤和正常样本的结果，筛选肿瘤特异性的 SV
```

## 移动元件插入（MEI）检测

```bash
# 使用 -q 参数检测 MEI
pindel-tool pindel -f reference.fa -i config.txt -o mei_output -c ALL -q

# 可选：使用 BreakDancer 结果辅助
pindel-tool pindel -f reference.fa -i config.txt -o mei_output -c ALL -q -b breakdancer_output.txt
```

## 常见问题

### 1. 内存使用过高

**原因**：分析着丝粒等高重复区域会消耗大量内存

**解决方案**：
```bash
# 排除问题区域
pindel-tool pindel -f ref.fa -i config.txt -o output -c ALL -J centromere.bed
```

### 2. 运行速度慢

**优化方案**：
- 使用多线程：`-T 4`
- 并行处理染色体
- 排除低价值区域：`-J exclude.bed`

### 3. BAM 文件格式错误

确保 BAM 文件：
- 已按参考基因组排序
- 已建立索引（`.bai` 文件）
- 包含正确的头信息

## ANNOVAR 注释

### 内含子变异 cDNA 坐标注释

Pindel 检测的结构变异可能位于内含子区域。本仓库修改了 ANNOVAR 源代码，使内含子变异的注释自动输出到 `AAChange.refGeneWithVer` 列，与外显子变异格式保持一致。

**特性**：
- 自动计算内含子变异的 cDNA 坐标
- 只保留每个转录本最近的外显子注释
- 输出格式添加基因名前缀

```bash
# 运行 ANNOVAR 注释（自动启用内含子 cDNA 坐标计算）
pindel-tool anno input.vcf /path/to/humandb/ \
    --buildver hg38 \
    --protocol refGeneWithVer \
    --nastring . \
    --out output
```

### 输出格式

| 列名 | 内含子变异 | 外显子变异 |
|------|-----------|-----------|
| Func.refGeneWithVer | intron | exonic |
| Gene.refGeneWithVer | 基因名 | 基因名 |
| GeneDetail.refGeneWithVer | . | cDNA 注释 |
| ExonicFunc.refGeneWithVer | . | 功能类型 |
| AAChange.refGeneWithVer | cDNA 注释 | 氨基酸改变 |

### 注释格式说明

内含子变异注释格式：`DEAF1:NM_021008.4:exon11:c.1504-2107->AAAAA`

| 字段 | 含义 |
|------|------|
| DEAF1 | 基因名 |
| NM_021008.4 | 转录本 ID |
| exon11 | 最近的外显子 |
| c.1504 | 外显子末端的 cDNA 位置 |
| -2107 | 内含子中的距离（负值=下游，正值=上游） |
| ->AAAAA | 变异类型（插入/删除/替换） |

### 源代码修改

修改了 `annovar_scripts/table_annovar.pl`：
1. 将内含子变异的 cDNA 注释从 `GeneDetail` 列移动到 `AAChange` 列
2. 添加 `filter_intronic_annotation` 函数，筛选每个转录本最近的外显子注释
3. 添加基因名前缀到注释输出

## 版本信息

- Pindel version: 0.2.5b9
- 发布日期: 20160729
- 作者: Kai Ye (kaiye@xjtu.edu.cn)
- Python CLI tool: pindel-tool 0.1.0（`pip install pindel-tool`）

## 参考文献

Ye K, Schulz MH, Long Q, Apweiler R, Ning Z. Pindel: a pattern growth approach to detect break points of large deletions and medium sized insertions from paired-end short reads. Bioinformatics. 2009;25(21):2865-2871.

## 许可证

详见 COPYING.txt 文件。

## 联系方式

- 问题反馈: https://github.com/genome/pindel/issues
- 主要作者: Kai Ye <kaiye@xjtu.edu.cn>
- pindel2vcf 相关: Eric-Wubbo Lameijer <e.m.w.lameijer@gmail.com>