# Pindel

Pindel 是一款基于 split-read（分裂读）算法的结构变异（Structural Variant, SV）检测工具，专门用于检测：

- **缺失 (Deletions, DEL)**
- **短插入 (Short Insertions, SI)**
- **长插入 (Long Insertions, LI)**
- **倒位 (Inversions, INV)**
- **串联重复 (Tandem Duplications, TD)**
- **移动元件插入 (Mobile Element Insertions, MEI)**

## 目录结构

```
pindel/
├── pindel              # 主程序
├── pindel2vcf          # VCF 转换工具
├── sam2pindel          # SAM 转换工具
├── pindel2vcf4tcga     # TCGA 格式转换工具
└── htslib/             # 本地 htslib 库
    ├── libhts.so.3
    └── htslib/         # 头文件
```

## 编译安装

### 前置要求

- GNU Make 和 GCC
- htslib（已包含在本目录的 `htslib/` 子目录中）

### 编译步骤

```bash
# 直接编译（使用本地 htslib）
make clean
make

# 编译产物
# - pindel: 主程序
# - pindel2vcf: VCF 转换工具
# - sam2pindel: SAM 转换工具
# - pindel2vcf4tcga: TCGA 格式转换工具
```

### GCC 13+ 兼容性

如果使用 GCC 13 或更高版本编译，源码已修复 `abs(unsigned int)` 类型歧义问题。

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
./pindel -f reference.fa -i config.txt -o output_prefix -c ALL

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
./pindel -f ref.fa -i config.txt -o output -c 20

# 分析染色体特定区域
./pindel -f ref.fa -i config.txt -o output -c 20:1000000-2000000

# 排除着丝粒/端粒区域（减少假阳性）
./pindel -f ref.fa -i config.txt -o output -c ALL -J centromere.bed
```

### 4. 性能优化

#### 多线程加速

```bash
# 使用 4 线程
./pindel -f ref.fa -i config.txt -o output -c ALL -T 4
```

#### 并行处理（集群环境）

```bash
# 按染色体并行提交任务
for chr in {1..22} X Y; do
    ./pindel -f ref.fa -i config.txt -o chr${chr} -c ${chr} &
done
wait
```

#### 降低内存使用

```bash
# 减小窗口大小
./pindel -f ref.fa -i config.txt -o output -c ALL -W 1

# 排除高重复区域
./pindel -f ref.fa -i config.txt -o output -c ALL -J exclude.bed
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

使用 `pindel2vcf` 将 Pindel 输出转换为 VCF：

```bash
./pindel2vcf -r reference.fa -R GRCh38 -d 20200101 -p output_D -e 5

# 参数说明：
# -r  参考基因组 FASTA 文件
# -R  参考基因组名称
# -d  参考基因组版本日期
# -p  Pindel 输出文件（如 output_D）
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
./pindel -f /data/reference.fa -i config.txt -o pindel_result -c ALL -T 4 -J centromere.bed

# 3. 转换为 VCF（过滤支持读数 <5 的结果）
./pindel2vcf -r /data/reference.fa -R GRCh38 -d 20200101 -p pindel_result_D -e 5 -v deletions.vcf
./pindel2vcf -r /data/reference.fa -R GRCh38 -d 20200101 -p pindel_result_INV -e 5 -v inversions.vcf
```

## 肿瘤-正常配对分析

对于肿瘤样本的体细胞突变检测：

```bash
# 创建配置文件
echo "/data/tumor.bam    350    TUMOR" > tn_config.txt
echo "/data/normal.bam   350    NORMAL" >> tn_config.txt

# 运行 Pindel
./pindel -f hg38.fa -i tn_config.txt -o somatic_sv -c ALL -T 4 -g

# 后续需要比较肿瘤和正常样本的结果，筛选肿瘤特异性的 SV
```

## 移动元件插入（MEI）检测

```bash
# 使用 -q 参数检测 MEI
./pindel -f reference.fa -i config.txt -o mei_output -c ALL -q

# 可选：使用 BreakDancer 结果辅助
./pindel -f reference.fa -i config.txt -o mei_output -c ALL -q -b breakdancer_output.txt
```

## 常见问题

### 1. 内存使用过高

**原因**：分析着丝粒等高重复区域会消耗大量内存

**解决方案**：
```bash
# 排除问题区域
./pindel -f ref.fa -i config.txt -o output -c ALL -J centromere.bed
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

## 版本信息

- Pindel version: 0.2.5b9
- 发布日期: 20160729
- 作者: Kai Ye (kaiye@xjtu.edu.cn)

## 参考文献

Ye K, Schulz MH, Long Q, Apweiler R, Ning Z. Pindel: a pattern growth approach to detect break points of large deletions and medium sized insertions from paired-end short reads. Bioinformatics. 2009;25(21):2865-2871.

## 许可证

详见 COPYING.txt 文件。

## 联系方式

- 问题反馈: https://github.com/genome/pindel/issues
- 主要作者: Kai Ye <kaiye@xjtu.edu.cn>
- pindel2vcf 相关: Eric-Wubbo Lameijer <e.m.w.lameijer@gmail.com>