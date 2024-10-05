from argparse   import ArgumentParser
from datetime   import datetime
from matplotlib import patches, pyplot, ticker
from numpy      import log10, max, min
from os         import mkdir, path
from pandas     import DataFrame, read_csv
from random     import choice

import re

gene_metrics        = ["50%-assembled genes",
                       "95%-assembled genes",
                       "50%-covered genes",
                       "95%-covered genes"]       
isoforms_metrics    = ["50%-assembled isoforms",
                        "95%-assembled isoforms",
                        "50%-covered isoforms",
                        "95%-covered isoforms"]
transcripts_metrics = ["Transcripts > 500 bp",
                        "Transcripts > 1000 bp",
                        "Aligned",
                        "Uniquely aligned",
                        "Multiply aligned",
                        "Unaligned",
                        "Misassemblies",
                        "Unannotated",
                        "50%-matched",
                        "95%-matched"]
scaled_metrics =       ["Mean fraction of transcript matched",
                        "Mean isoform assembly",
                        "Mean isoform coverage",
                        "Database coverage",
                        "Avg. aligned fraction",
                        ]

def scientific_format(x, pos):
        if x == 0:
            return '0'
        exponent = int(log10(x))
        base = x / (10**exponent)
        return f'{base:.0f}*10^{exponent}'

class MyArgumentParser(ArgumentParser):

    prog        =   "rnaQUASTplotter"

    description =   """
                    ----------\n
                    Commandline tool for comparative plots for rnaQUAST reports
                    for multiple assemblies\n
                    ----------
                    """
    
    help = {
            "reportdirs": "paths to output directories from rnaQUAST",
            "names":      "list of names for the assemblies (default=[\"auto\"])",
            "colors":     "list of colors in hexcode (default=[\"auto\"])",
        }
    
    def __init__(self) -> None:

        super().__init__(prog=self.prog, description=self.description)
        
        self.add_argument("report_dirs", nargs='+', type=str, help=self.help["reportdirs"])
        self.add_argument("-names",      nargs='+', type=str, help=self.help["names"],      default=["auto"])
        self.add_argument("-colors",     nargs='+', type=str, help=self.help["colors"],     default=["auto"])

class ReportParser():

    @classmethod
    def find_number(cls, string: str)->float:
        """
        Parameters:
            string (str)
        Returns:
            float, the first number found within the string

        A method used to extract the values from the database_metrics.txt-file
        """
        # Regular expression to find the first number (integer or float) in the string
        match = re.search(r'[-+]?\d*\.?\d+', string)
        
        return float(match.group(0))

    @classmethod
    def parse_database_metrics(cls, file_path: str)->dict:
        """
        Parameters:
            file_path (str), path of a database_metrics.txt-file
        Returns:
            dict, a dictionary with the number of genes and isoforms in the database
        """
        with open(file_path, "r") as file:
            lines = file.readlines()

        return {"Genes": cls.find_number(lines[1]), "Isoforms": cls.find_number(lines[4])}

    @classmethod
    def parse_report(cls, report_dir):

        short_report = read_csv(path.join(report_dir, "short_report.tsv"), sep="\t")
        
        short_report.rename(columns={'METRICS/TRANSCRIPTS': 'metrics'}, inplace=True)

        file_path = path.join(report_dir, short_report.columns[1]+"_output")
        file_path = path.join(file_path, "database_metrics.txt")

        database_metrics = cls.parse_database_metrics(file_path)

        return short_report, database_metrics
    
    @classmethod
    def parse_reports(cls, report_dirs):

        short_reports    = []
        database_metrics = []

        for report_dir in report_dirs:

            sr, dm = cls.parse_report(report_dir)
            short_reports.append(sr)
            database_metrics.append(dm)

        if len([dm for dm in database_metrics if not database_metrics[0]["Genes"]==dm["Genes"]]):
            raise Exception("Database metrics vary for some assemblies")
        else:
            return short_reports, database_metrics[0]
        
class ValueScaler():
    
    @classmethod
    def mmpt_to_mmpkb(cls, short_report: DataFrame):
        """
        Converts the "Avg. mismatches per transcript" to "Avg. mismatches per aligned kb"
        """
        sr = short_report
        avg_alignment_len = float(sr.loc[sr["metrics"]=="Avg. alignment length"][sr.columns[1]])
        avg_mismatches    = float(sr.loc[sr["metrics"]=="Avg. mismatches per transcript"][sr.columns[1]])
        i = sr.loc[sr["metrics"]=="Avg. mismatches per transcript"].index
        sr.loc[i,"metrics"]     = "Avg. mismatches per aligned kb"
        sr.loc[i,sr.columns[1]] = 1000 * avg_mismatches / avg_alignment_len
    
    @classmethod
    def find_divider(cls, short_reports: list, i: int, metric: str, database_metrics: DataFrame):

        if metric in gene_metrics:
            return database_metrics["Genes"]
        if metric in isoforms_metrics:
            return database_metrics["Isoforms"] 
        if metric in transcripts_metrics:
            sr = short_reports[i]
            return sr.loc[sr["metrics"]=="Transcripts"][sr.columns[1]]
        if metric in scaled_metrics:
            return 1
        return max([float(sr.loc[sr["metrics"]==metric][sr.columns[1]]) for sr in short_reports])
    
    @classmethod
    def scale(cls, short_reports, database_metrics):

        scaled_values = [[] for sr in short_reports]

        for i, sr in enumerate(short_reports):
            cls.mmpt_to_mmpkb(sr)
        
        for i, sr in enumerate(short_reports):
            for j, metric in enumerate(sr["metrics"]):
                value   = sr.loc[sr["metrics"]==metric][sr.columns[1]]
                divider = cls.find_divider(short_reports, i, metric, database_metrics)
                scaled_values[i].append(float(value) / float(divider))
            
        for i, sr in enumerate(short_reports):
            sr['scaled'] = scaled_values[i]
        
class Plotter():

    @classmethod
    def empty_plot(cls, axes, metrics, names, colors):

        for tick in axes.get_xticks():
            axes.axvline(x=tick, color='gray', linestyle='--', alpha=0.5)
        axes.tick_params(axis='x', which='major', labelsize=30)
        axes.set_yticks(list(range(len(metrics))))
        axes.set_yticklabels(metrics, fontweight='bold', fontsize=30)
        axes.set_ylim(ymin=-0.5,ymax=len(metrics)-0.5)

    @classmethod
    def add_legend(cls, axes, names, colors):

        handles = [patches.Rectangle([0,0],5,5,color=c) for c in colors]
        axes.legend(handles, names, bbox_to_anchor=(-2, -0.9, 1, 1), fontsize=30)

    @classmethod
    def fill_plot_bars(cls, axes, short_reports, metrics, colors, scaled=False):

        width  = 1 / len(short_reports)
        offset = width / 2 - 0.5
        xmin   = 1
        xmax   = 0
        
        for i, metric in enumerate(metrics):
            for j, sr in enumerate(short_reports):

                if scaled:
                    col = "scaled"
                else:
                    col = sr.columns[1]

                y     = offset + i + (j * width)
                x     = sr.loc[sr["metrics"]==metric][col]
                xmax = max([float(x), xmax])
                xmin = min([float(x),xmin])
                axes.barh(y, x, color=colors[j], height=width)

        axes.set_xlim(xmin=0,xmax=xmax)
        for tick in axes.get_xticks():
            axes.axvline(x=tick, color='gray', linestyle='--', alpha=0.5)

        axes.set_xticklabels(axes.get_xticks(), rotation=45)

    @classmethod
    def fill_plot_lines(cls, axes, short_reports, metrics, colors, scaled=False):

        xmin = 1
        xmax = 0

        for i,sr in enumerate(short_reports):

            if scaled:
                col = "scaled"
            else:
                col = sr.columns[1]

            x = sr.loc[sr["metrics"].isin(metrics)][col]
            xmax = max([max(x),xmax])
            xmin = min([min(x),xmin])
            axes.plot(x,list(range(len(x))), c=colors[i], linewidth=5)
        
        axes.set_xlim(xmin=0,xmax=xmax)
        for tick in axes.get_xticks():
            axes.axvline(x=tick, color='gray', linestyle='--', alpha=0.5)

        axes.set_xticklabels(axes.get_xticks(), rotation=45)

    @classmethod
    def generate_plots(cls, short_reports, names, colors, save_as, n_isoforms):

        # 1 Collective scaled plots

        metrics = short_reports[0]['metrics'][2:]

        # 1.1 Line plot 
        fig, axes = pyplot.subplots(figsize=(5,len(metrics)))
        cls.empty_plot(axes, metrics, names, colors)
        cls.fill_plot_lines(axes, short_reports, metrics, colors, scaled=True)
        cls.add_legend(axes, names, colors)
        pyplot.savefig(save_as+"_scaled_lines", bbox_inches='tight', pad_inches=0.5)

        # 1.2 Bar plot
        fig, axes = pyplot.subplots(figsize=(5,len(metrics)))
        cls.empty_plot(axes, metrics, names, colors)
        cls.fill_plot_bars(axes, short_reports, metrics, colors, scaled=True)
        cls.add_legend(axes, names, colors)
        pyplot.savefig(save_as+"_scaled_bars", bbox_inches='tight', pad_inches=0.5)

        # 2. Unscaled group plots

        for sr in short_reports:
            sr.loc[len(sr)] = ["Isoforms", n_isoforms, 1]

        metrics_list = [["Genes"]+gene_metrics,
                        ["Transcripts"]+transcripts_metrics,
                        ["Isoforms"]+isoforms_metrics,
                        ["Avg. aligned fraction",
                         "Avg. mismatches per aligned kb",
                         "Database coverage",
                         "Duplication ratio",
                         "Mean isoform coverage",
                         "Mean isoform assembly",
                         "Mean fraction of transcript matched"]]
        group_names  = ["Gene metrics",
                        "Transcript metrics",
                        "Isoform metrics",
                        "Other metrics"]
        
        n = len(group_names)

        for i in range(n):

            fig, axes = pyplot.subplots(figsize=(5,len(metrics)/2))

            cls.empty_plot(axes, metrics_list[i], names, colors)
            cls.fill_plot_lines(axes, short_reports, metrics_list[i], colors, scaled=False)
            axes.set_title(group_names[i], fontweight="bold", fontsize=35)
        
            pyplot.savefig(save_as+f"_absolute_lines_{group_names[i]}_no_legend", bbox_inches='tight', pad_inches=0.5)
            cls.add_legend(axes, names, colors)
            pyplot.savefig(save_as+f"_absolute_lines_{group_names[i]}_legend", bbox_inches='tight', pad_inches=0.5)

        for i in range(n):

            fig, axes = pyplot.subplots(figsize=(5,len(metrics)/2))

            cls.empty_plot(axes, metrics_list[i], names, colors)
            cls.fill_plot_bars(axes, short_reports, metrics_list[i], colors, scaled=False)
            axes.set_title(group_names[i], fontweight="bold", fontsize=35)

            pyplot.savefig(save_as+f"_absolute_bars_{group_names[i]}_no_legend", bbox_inches='tight', pad_inches=0.5)
            cls.add_legend(axes, names, colors)
            pyplot.savefig(save_as+f"_absolute_bars_{group_names[i]}_legend", bbox_inches='tight', pad_inches=0.5)
            
def random_color():
        return "#" + "".join([choice("0123456789abcdef") for _ in range(6)])

def main():

    args = MyArgumentParser().parse_args()

    short_reports, database_metrics = ReportParser().parse_reports(args.report_dirs)

    # Sanity checks and automatically setting nonprovided optional parameters
    if args.names == ["auto"]:
        names  = [sr.columns[1] for sr in short_reports]
    else:
        if not len(args.report_dirs) == len(args.names):
            raise Exception("Number of names must match number of reports")
        else:
            names = args.names
    if args.colors == ["auto"]:
        colors = [random_color() for _ in short_reports]
    else:
        if not len(args.report_dirs) == len(args.names):
            raise Exception("Number of colors must match number of reports")
        else:
            colors = args.colors

    # Data processing
    ValueScaler().scale(short_reports, database_metrics)

    # Prepare output
    out_dir = datetime.now().strftime("%dd%mm%Yy_%Hh%Mm%Ss")
    save_as = path.join(out_dir, "rnaQUAST_comparison")
    mkdir(out_dir)

    # Combined dataframe
    combined = DataFrame({"metrics": short_reports[0]['metrics']})
    for i, sr in enumerate(short_reports):
        combined[names[i]+' (absolute)'] = sr[sr.columns[1]]
        combined[names[i]+' (scaled)']   = sr[sr.columns[2]]  

    print(combined)

    combined.to_csv(save_as+".csv")
    combined.to_csv(save_as+".tsv", sep="\t")

    for i in range(combined.shape[0]):
        combined.loc[i,"metrics"] = str(combined.loc[i,"metrics"]).replace("%","\\%")

    combined.style.to_latex(save_as+".tex")

    # Plotting
    Plotter().generate_plots(short_reports, names, colors, save_as, database_metrics["Isoforms"])

    print("Done")
    exit(0)

if __name__ == '__main__':
    
    main()
        