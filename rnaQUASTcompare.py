from argparse   import ArgumentParser
from datetime   import datetime
from matplotlib import patches, pyplot
from numpy      import max
from os         import mkdir, path
from pandas     import DataFrame, read_csv
from random     import choice

import re

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

    gene_metrics        = ["50%-assembled genes",
                            "95%-assembled genes",
                            "50%-covered genes",
                            "95%-covered genes"]       
    isoforms_metrics    = ["50%-assembled isoforms",
                            "95%-assembled isoforms",
                            "50%-covered isoforms",
                            "50%-covered isoforms"]
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
    scaled_metrics =       ["Database coverage",
                            "Avg. aligned fraction",
                            "Mean fraction of transcript matched"]
    
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

        if metric in cls.gene_metrics:
            return database_metrics["Genes"]
        if metric in cls.isoforms_metrics:
            return database_metrics["Isoforms"] 
        if metric in cls.transcripts_metrics:
            sr = short_reports[i]
            return sr.loc[sr["metrics"]=="Transcripts"][sr.columns[1]]
        if metric in cls.scaled_metrics:
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
    def empty_plot(cls, metrics, names, colors):

        fig, axes = pyplot.subplots(figsize=(5,5*len(names)))

        for tick in axes.get_xticks():
            axes.axvline(x=tick, color='gray', linestyle='--', alpha=0.5)
        axes.tick_params(axis='x', which='major', labelsize=20)
        axes.set_yticks(list(range(len(metrics)-2)))
        axes.set_yticklabels(metrics[2:], fontweight='bold', fontsize=20)
        axes.set_xlim(xmin=0,xmax=1)
        axes.set_ylim(ymin=-0.5,ymax=len(metrics)-2.5)

        handles = [patches.Rectangle([0,0],5,5,color=c) for c in colors]

        axes.legend(handles, names, bbox_to_anchor=(0.5, -0.55, 0.5, 0.5), fontsize=20)

        return fig, axes

    @classmethod
    def fill_plot_bars(cls, axes, short_reports, colors):

        width  = 1 / len(short_reports)
        offset = width / 2 - 0.5
        
        for i, metric in enumerate(short_reports[0]['metrics'][2:]):

            for j, sr in enumerate(short_reports):

                y = offset + i + (j * width)
                x = sr.loc[sr["metrics"]==metric]["scaled"]
                axes.barh(y, x, color=colors[j], height=width)

    @classmethod
    def fill_plot_lines(cls, axes, short_reports, colors):

        for i,sr in enumerate(short_reports):
            axes.plot(sr['scaled'][2:],list(range(sr.shape[0]-2)), c=colors[i], linewidth=5)

    @classmethod
    def plot(cls, short_reports, names, colors, save_as):

        # Line plot
        fig, axes = cls.empty_plot(short_reports[0]["metrics"], names, colors)
        cls.fill_plot_lines(axes, short_reports, colors)
        pyplot.savefig(save_as+"_lines", bbox_inches='tight', pad_inches=0.5)

        # Bar plot
        fig, axes = cls.empty_plot(short_reports[0]["metrics"], names, colors)
        cls.fill_plot_bars(axes, short_reports, colors)
        pyplot.savefig(save_as+"_bars", bbox_inches='tight', pad_inches=0.5)

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

    # Plotting
    Plotter().plot(short_reports, names, colors, save_as)

    # Combined dataframe
    combined = DataFrame({"metric": short_reports[0]['metrics']})
    for i, sr in enumerate(short_reports):
        combined[names[i]+' (absolute)'] = sr[sr.columns[1]]
        combined[names[i]+' (scaled)']   = sr[sr.columns[2]]
    print(combined)
    combined.to_csv(save_as+".tsv", sep="\t")
    for i in range(combined.shape[0]):
        combined.loc[i,"metric"] = str(combined.loc[i,"metric"]).replace("%","\\%")
    combined.style.to_latex(save_as+".tex")
    
    print("Done")
    exit(0)

if __name__ == '__main__':
    
    main()
        

        

