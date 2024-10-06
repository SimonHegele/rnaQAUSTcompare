# rnaQUASTcompare

Small python command line tool that generates a comparative plots for multiple rnaQUAST short reports.<br>
rnaQUAST (https://github.com/ablab/rnaquast) is a great tool for the evaluation of transcriptome assemblies.<br>

It generates a multitude of metrics for the quality of transcriptome assemblies, many of them by mapping the transcripts to an annotated genome.<br>
rnaQUAST does a great job at rating individual assemblies, however, directly comparing different reports is not as easy.

This tools compares metrics from rnaQUASTs short reports.

## Usage

```
positional arguments:
  report_dirs           paths to output directories from rnaQUAST

options:
  -h, --help            show this help message and exit
  -names NAMES [NAMES ...]
                        list of names for the assemblies (default=["auto"])
  -colors COLORS [COLORS ...]
                        list of colors in hexcode (default=["auto"])
```

## Output

rnaQUASTcompare.py will generate a folder with the current date and time in the same directory.<br>

I found the metric "Avg. mismatches per transcripts" to favor assemblies with transcripts that are shorter and replaced it with "Avg. mismatches per aligned kb".

### 1. Dataframes

Dataframes combining the data of all short reports in .csv, .tsv and .tex format

### 2. Plots

Metrics are grouped into four groups: "Gene metrics", "Transcript metrics", "Isoform metrics"<br>
and other metrics. For each of them a bar and a line plot will be created and a bar and line<br>
for all metrics together is created. In the comined plot all values are scaled to [0,1].<br>

Combined plots for all metrics with scaled values and individual plots for each metrics group.

**Example:**

A comparison of three transcriptome assembly tools from the same RNA-Seq data.

<p float="left">
  <img src="output_example/rnaQUAST_comparison_absolute_lines_Gene metrics_no_legend.png" width="400" />
  <img src="output_example/rnaQUAST_comparison_absolute_bars_Isoform metrics_no_legend.png" width="400" /> 
</p>

<p float="left">
  <img src="output_example/rnaQUAST_comparison_absolute_lines_Transcript metrics_legend.png" width="400" />
  <img src="output_example/rnaQUAST_comparison_absolute_bars_Other metrics_legend.png" width="400" /> 
</p>

**Value scaling**

I divided the metrics into groups:

- Gene metrics<br>
"50%-assembled genes", "95%-assembled genes", "50%-covered genes", "95%-covered genes"
- Isoforms metrics<br>
"50%-assembled isoforms", "95%-assembled isoforms", "50%-covered isoforms", "50%-covered isoforms"
- Transcript metrics<br>
"Transcripts > 500 bp", "Transcripts > 1000 bp", "Aligned", "Uniquely aligned", "Multiply aligned", "Unaligned", "Misassemblies", "Unannotated", "50%-matched", "95%-matched"
- Scaled metrics<br>
"Database coverage", "Avg. aligned fraction", "Mean fraction of transcript matched"
- Other metrics<br>
"Transcripts", "Avg. mismatches per aligned kb", "Duplication ratio"

Gene metrics are divided by the number of genes in the genome annotation<br>
Isoforms metrics are divided by the number of isoforms in the genome annotation.<br>
Isoforms metrics are divided by the number of sequences in the respective assembly.<br>
Scaled metrics are left unchanged.<br>
Other metrics are divided by the maximum value for all assemblies.
