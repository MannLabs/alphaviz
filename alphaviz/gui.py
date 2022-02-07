import os
import logging
import platform
import json
import pandas as pd
import numpy as np
from io import StringIO

import alphatims.bruker
import alphatims.utils

# visualization
import panel as pn
import bokeh.server.views.ws
from bokeh.models.widgets.tables import NumberFormatter
import plotly.express as px
import holoviews as hv
import matplotlib.colors

# local
import alphaviz
import alphaviz.utils
import alphaviz.io
import alphaviz.preprocessing
import alphaviz.plotting


def get_css_style(
    file_name="dashboard_style.css",
    directory=alphaviz.utils.STYLE_PATH
):
    file = os.path.join(
        directory,
        file_name
    )
    with open(file) as f:
        return f.read()


def init_panel():
    pn.extension(raw_css=[get_css_style()])
    hv.extension('bokeh')
    pn.extension('plotly')
    pn.extension('tabulator')


def update_config(filename, height=500, width=1500, ext='svg'):
    config = {
        'displaylogo': False,
        'toImageButtonOptions': {
            'format': f'{ext}', # one of png, svg, jpeg, webp
            'filename': f'{filename}',
            'height': height,
            'width': width,
            'scale': 1 # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    return config

if platform.system() == 'Windows':
    raw_folder_placeholder = r'D:\bruker\21min_HELA_proteomics'
    output_folder_placeholder = r'D:\bruker\21min_HELA_proteomics\txt'
    fasta_path_placeholder = r'D:\fasta_files\human.fasta'
else:
    # TODO: add a linux support
    raw_folder_placeholder = '/Users/test/bruker/21min_HELA_proteomics'
    output_folder_placeholder = '/Users/test/bruker/21min_HELA_proteomics/txt'
    fasta_path_placeholder = '/Users/test/fasta_files/human.fasta'


class BaseWidget(object):

    def __init__(self, name):
        self.name = name
        self.__update_event = pn.widgets.IntInput(value=0)
        self.depends = pn.depends(self.__update_event.param.value)
        self.active_depends = pn.depends(
            self.__update_event.param.value,
            watch=True
        )

    def trigger_dependancy(self):
        self.__update_event.value += 1


class HeaderWidget(object):
    """This class creates a layout for the header of the dashboard with the name of the tool and all links to the MPI website, the MPI Mann Lab page and the GitHub repo.

    Parameters
    ----------
    title : str
        The name of the tool.

    Attributes
    ----------
    header_title : pn.pane.Markdown
        A Panel Markdown pane that returns the title of the tool.
    mpi_biochem_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the MPI Biochemisty logo and makes the image clickable with the link to the official website.
    mpi_logo : pn.pane.JPG
        A Panel JPG pane that embeds a jpg image file of the MPI Biochemisty logo and makes the image clickable with the link to the official website.
    github_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the GitHub logo and makes the image clickable with the link to the GitHub repository of the project.

    """

    def __init__(
        self,
        title,
        img_folder_path,
        github_url
    ):
        self.layout = None
        self.header_title = pn.pane.Markdown(
            f'# {title}',
            sizing_mode='stretch_width',
        )
        self.biochem_logo_path = os.path.join(
            img_folder_path,
            "mpi_logo.png"
        )
        self.mpi_logo_path = os.path.join(
            img_folder_path,
            "max-planck-gesellschaft.jpg"
        )
        self.github_logo_path = os.path.join(
            img_folder_path,
            "github.png"
        )
        self.mpi_biochem_logo = pn.pane.PNG(
            self.biochem_logo_path,
            link_url='https://www.biochem.mpg.de/mann',
            width=60,
            height=60,
            align='start'
        )
        self.mpi_logo = pn.pane.JPG(
            self.mpi_logo_path,
            link_url='https://www.biochem.mpg.de/en',
            height=62,
            embed=True,
            width=62,
            margin=(5, 0, 0, 5),
            css_classes=['opt']
        )
        self.github_logo = pn.pane.PNG(
            self.github_logo_path,
            link_url=github_url,
            height=70,
            align='end'
        )

    def create_layout(self):
        self.layout = pn.Row(
            self.mpi_biochem_logo,
            self.mpi_logo,
            self.header_title,
            self.github_logo,
            height=73,
            sizing_mode='stretch_width'
        )
        return self.layout


class MainWidget(object):
    """This class create a layout for the main part of the dashboard with the description of the tool and a button to download the manual for the project's GUI.

    Parameters
    ----------
    description : str
        The short description of the tool.
    manual_path : str
        The path to the GUI manual.

    Attributes
    ----------
    project_description : pn.pane.Markdown
        A Panel Markdown pane that shows the description of the project.
    manual : pn.widgets.FileDownload
        A Panel FileDownload widget that allows to download the GUI manual of the tool.

    """
    def __init__(
        self,
        description,
        manual_path
    ):
        self.layout = None
        self.project_description = pn.pane.Markdown(
            description,
            margin=(10, 0, 10, 0),
            css_classes=['main-part'],
            align='start',
            width=460
        )
        self.manual = pn.widgets.FileDownload(
            file=manual_path,
            label='Download Manual',
            button_type='default',
            align='center',
            auto=True,
            height=31,
            width=200,
            margin=(0, 20, 0, 0)
        )
        self.download_new_version_button = pn.widgets.Button(
            button_type="danger",
            align='center',
            height=31,
            width=200,
            margin=(20, 20, 0, 0)
        )

    def create_layout(self):

        latest_github_version = alphaviz.utils.check_github_version(silent=False)

        if latest_github_version and latest_github_version != alphaviz.__version__:
            self.download_new_version_button.name = f"Download version {latest_github_version}"
            download_new_version_button = self.download_new_version_button
            download_new_version_button.js_on_click(
                code=f"""window.open("https://github.com/MannLabs/alphaviz/releases/latest")"""
            )
        else:
            download_new_version_button = None

        self.layout = pn.Row(
            self.project_description,
            pn.layout.HSpacer(width=500),
            pn.Column(
                self.manual,
                download_new_version_button,
                align='center',
            ),
            background='#eaeaea',
            align='center',
            sizing_mode='stretch_width',
            height=190,
            margin=(10, 8, 10, 8),
            css_classes=['background']
        )
        return self.layout


class DataImportWidget(BaseWidget):

    def __init__(self):
        super().__init__(name="Data")
        self.raw_data = None
        self.mq_evidence = None
        self.mq_all_peptides = None
        self.mq_msms = None
        self.mq_protein_groups = None
        self.mq_summary = None
        self.diann_proteins = None
        self.diann_peptides = None
        self.diann_statist = None
        self.fasta = None
        self.layout = None
        self.settings = {
            'path_evidence_file': str(),
            'analysis_software': str()
        }
        self.path_raw_folder = pn.widgets.TextInput(
            name='Specify the full path to the folder with unprocessed Bruker files:\u002a',
            placeholder=raw_folder_placeholder,
            width=900,
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.ms_file_name = pn.widgets.Select(
            name='Select the raw file:\u002a',
            size=10,
            width=900,
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.path_output_folder = pn.widgets.TextInput(
            name='Specify the full path to the output folder of any supported software analysis tool:',
            placeholder=output_folder_placeholder,
            width=900,
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.path_fasta_file = pn.widgets.TextInput(
            # TODO: remove the fixed fasta file before release
            # value='/Users/eugeniavoytik/copied/Bruker/MaxQuant_output_tables/20210413_TIMS03_EVO03_PaSk_MA_HeLa_200ng_S1-A1_1_24848.d/txt/human.fasta',
            name='Specify the full path to the fasta file:',
            placeholder=fasta_path_placeholder,
            width=900,
            sizing_mode='stretch_width',
            margin=(15, 15, 15, 15)
        )
        self.is_prediction = pn.widgets.Checkbox(
            name='Activate the prediction',
            margin=(5, 0, 5, 15),
            disabled=True,
        )
        # UPLOAD DATA
        self.upload_button = pn.widgets.Button(
            name='Load Data',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.upload_progress = pn.indicators.Progress(
            max=1,
            value=1,
            active=False,
            bar_color='light',
            width=250,
            align='center',
            margin=(-10, 0, 30, 0)
        )
        self.import_error = pn.pane.Alert(
            alert_type="danger",
            object='',
            sizing_mode='stretch_width',
            margin=(10, 0, 5, 0),
        )

    def create_layout(self):
        dependances = {
            self.path_raw_folder: [self.update_file_names, 'value'],
            self.ms_file_name: [self.update_output_folder_and_fasta, 'value'],
            self.upload_button: [self.load_data, 'clicks'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1]
            )
        self.layout = pn.Card(
            pn.Row(
                pn.Column(
                    self.path_raw_folder,
                    self.ms_file_name,
                    self.path_output_folder,
                    self.path_fasta_file,
                    self.is_prediction,
                    margin=(10, 30, 10, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                pn.Column(
                    self.upload_button,
                    self.upload_progress,
                    self.import_error,
                    align='center',
                    margin=(100, 40, 0, 0),
                )
            ),
            title='Data Import',
            collapsed=False,
            collapsible=True,
            header_background='#eaeaea',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )
        return self.layout

    def update_file_names(self, *args):
        try:
            self.ms_file_name.options = alphaviz.preprocessing.sort_naturally(
                alphaviz.io.get_filenames_from_directory(
                    self.path_raw_folder.value,
                    ['d', 'hdf']
                )
            )
        except OSError:
            self.import_error.object = f"#### The selected directory is not found."

    def update_output_folder_and_fasta(self, *args):
        try:
            data_folder = os.path.join(
                self.path_raw_folder.value,
                self.ms_file_name.value
            )
            for dirpath, dirnames, filenames in os.walk(data_folder):
                for filename in filenames:
                    if filename.endswith(".fasta"):
                        self.path_fasta_file.value = os.path.join(dirpath, filename)
                    elif filename == 'evidence.txt':
                        self.path_output_folder.value = dirpath
            if not self.path_fasta_file.value:
                for filename in os.listdir(self.path_raw_folder.value):
                    if filename.endswith(".fasta"):
                        self.path_fasta_file.value = os.path.join(self.path_raw_folder.value, filename)
        except:
            self.import_error.object = "#### The selected folder does not contain any .d or .hdf files."

    def load_data(self, *args):
        alphatims.utils.set_progress_callback(self.upload_progress)
        self.settings['analysis_software'] = ''
        self.import_error.object = ''
        self.upload_progress.value = 0
        try:
            self.raw_data = alphatims.bruker.TimsTOF(
                os.path.join(
                    self.path_raw_folder.value,
                    self.ms_file_name.value
                ),
            )
        except:
            self.import_error.object += '\n#### The selected unprocessed Bruker file is corrupted and cannot be loaded. \n#### Please select another file.',
            raise OSError('The selected unprocessed Bruker file is corrupted and cannot be loaded. Please select another file.')
        alphatims.utils.set_progress_callback(True)
        # TODO: to change it just by changing the active=True parameter for the self.upload_progress when the bug will be fixed
        self.upload_progress = pn.indicators.Progress(
            active=True,
            bar_color='light',
            width=250,
            align='center',
            margin=(-10, 0, 30, 0)
        )
        self.layout[0][2][1] = self.upload_progress

        # read the fasta file if specified
        if self.path_fasta_file.value:
            try:
                self.fasta = alphaviz.io.read_fasta(
                    self.path_fasta_file.value
                )
            except:
                self.import_error.object += "\n#### The selected fasta file cannot be loaded."
        else:
            self.import_error.object += "\n#### The fasta file file has not been provided."

        # read analysis output files (MQ, DIA-NN, etc.) if specified
        ## check all files in the analysis output folder
        if self.path_output_folder.value:
            files = alphaviz.io.get_filenames_from_directory(
                directory=self.path_output_folder.value,
                extensions_list=['txt', 'tsv', 'csv']
            )
            mq_files = ['allPeptides.txt', 'msms.txt', 'evidence.txt', 'proteinGroups.txt', 'summary.txt']
            if any(file in files for file in mq_files):
                print('Reading the MaxQuant output files...')
                if all(file in files for file in mq_files):
                    self.mq_all_peptides, self.mq_msms, self.mq_evidence, self.mq_protein_groups, self.mq_summary = alphaviz.io.import_mq_output(
                        mq_files,
                        self.path_output_folder.value,
                        self.ms_file_name.value.split('.')[0]
                    )
                    self.settings['analysis_software'] = 'maxquant'
                else:
                    self.import_error.object += "\n#### The MQ output files necessary for the visualization are not found."
            else:
                print('Reading the DIA-NN output files...')
                try:
                    self.diann_proteins, self.diann_peptides, self.diann_statist = alphaviz.io.import_diann_output(
                        self.path_output_folder.value,
                        self.ms_file_name.value.split('.')[0],
                        self.fasta
                    )
                    self.settings['analysis_software'] = 'diann'
                except:
                    self.import_error.object += "\n#### The DIA-NN output files necessary for the visualization are not found."
        else:
            self.import_error.object += "\n#### The output files of the supported software tools have not been provided."

        self.trigger_dependancy()
        self.upload_progress.active = False
        self.upload_progress.value = 100


class OptionsWidget(object):

    def __init__(self, data):
        self.data = data
        self.layout = pn.Card(
            title='Settings',
            collapsed=True,
            sizing_mode='stretch_width',
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

    def get_layout(self, *args):
        return self.layout

    def add_option(self, option):
        self.layout.append(option)


class HeatmapOptionsWidget(object):

    def __init__(self):
        self.plot1_x_axis = pn.widgets.Select(
            name='X-axis label',
            value='m/z, Th',
            options=['m/z, Th', 'Inversed IM, V·s·cm\u207B\u00B2'],
            width=180,
            margin=(20, 20, 20, 20),
        )
        self.plot1_y_axis = pn.widgets.Select(
            name='Y-axis label',
            value='Inversed IM, V·s·cm\u207B\u00B2',
            options=['m/z, Th', 'Inversed IM, V·s·cm\u207B\u00B2'],
            width=180,
            margin=(20, 20, 20, 10),
        )
        self.heatmap_color = pn.widgets.Select(
            name='Color scale',
            value='fire',
            options=sorted(hv.plotting.util.list_cmaps(reverse=True) + hv.plotting.util.list_cmaps(reverse=False)),
            width=180,
            margin=(20, 20, 20, 10),
        )
        self.heatmap_background = pn.widgets.Select(
            name='Background color',
            value='black',
            options=list(matplotlib.colors.CSS4_COLORS.keys()),
            width=180,
            margin=(20, 20, 20, 10),
        )
        self.precursor_target_size = pn.widgets.IntInput(
            name='Precursor target size',
            value=20,
            step=5,
            start=0,
            end=100,
            width=180,
            margin=(20, 20, 20, 10),
        )
        self.precursor_target_color = pn.widgets.Select(
            name='Precursor target color',
            value='lime',
            options=list(matplotlib.colors.CSS4_COLORS.keys()),
            width=180,
            margin=(20, 20, 20, 10),
        )

    def create_layout(self, *args):
        layout = pn.Card(
            pn.Row(
                self.plot1_x_axis,
                self.plot1_y_axis,
                self.heatmap_color,
                self.heatmap_background,
                self.precursor_target_size,
                self.precursor_target_color
            ),
            title='Heatmap options',
            collapsed=False,
            sizing_mode='stretch_width',
            margin=(15, 8, 0, 8),
            css_classes=['background']
        )
        return layout


class ToleranceOptionsWidget(object):

    def __init__(self):
        self.mz_tolerance = pn.widgets.FloatInput(
            name='m/z Tolerance (ppm)',
            value=10,
            step=5,
            start=0,
            end=1000,
            width=150,
            margin=(20, 20, 20, 10),
        )
        self.im_tolerance = pn.widgets.FloatInput(
            name='IM Tolerance (1/K0)',
            value=0.05,
            step=0.1,
            start=0,
            end=2,
            width=150,
            margin=(20, 20, 20, 10),
        )
        self.rt_tolerance = pn.widgets.FloatInput(
            name='RT Tolerance (sec)',
            value=30,
            step=5,
            start=0,
            end=1000,
            width=150,
            margin=(20, 20, 20, 10),
        )

    def create_layout(self, *args):
        layout = pn.Card(
            pn.Row(
                self.mz_tolerance,
                self.im_tolerance,
                self.rt_tolerance
            ),
            title='Tolerance settings',
            collapsed=False,
            sizing_mode='stretch_width',
            margin=(15, 8, 0, 8),
            css_classes=['background']
        )
        return layout

class ColorscaleOptionsWidget(object):

    def __init__(self):
        self.colorscale_qualitative = pn.widgets.Select(
            name='Qualitative color scale',
            value='Alphabet',
            options=list(set([each['y'][0] for each in px.colors.qualitative.swatches()['data']])),
            width=180,
            margin=(20, 20, 20, 10),
        )
        self.colorscale_sequential = pn.widgets.Select(
            name='Sequential color scale',
            value='Viridis',
            options=list(set([each['y'][0] for each in px.colors.sequential.swatches()['data']])),
            width=190,
            margin=(20, 20, 20, 10),
        )

    def create_layout(self, *args):
        layout = pn.Card(
            pn.Row(
                self.colorscale_qualitative,
                self.colorscale_sequential,
            ),
            title='Customization options',
            collapsed=False,
            sizing_mode='stretch_width',
            margin=(15, 8, 15, 8),
            css_classes=['background']
        )
        return layout


class TabsWidget(object):

    def __init__(self, data, options):
        self.layout = None
        self.data = data
        self.options = options

    def create_layout(
        self,
        tab_list=None
    ):
        self.tabs = tab_list
        return self.data.depends(self.return_layout)

    def return_layout(self, *args):
        if self.data.raw_data is not None:
            self.layout = pn.Tabs(
                tabs_location='above',
                margin=(10, 10, 5, 8),
                sizing_mode='stretch_width',
            )
            self.layout += self.tabs
            self.layout[0] = (
                'Main View',
                MainTab(self.data, self.options).create_layout()
            )
            self.layout[1] = (
                'Quality Control',
                QCTab(self.data, self.options).create_layout()
            )
            self.layout[2] = (
                'Targeted Mode',
                TargetModeTab(self.data, self.options).create_layout()
            )
            self.active = 0
            self.data.layout.collapsed = True
        return self.layout


class MainTab(object):

    def __init__(self, data, options):
        self.data = data
        self.mass_dict = alphaviz.utils.get_mass_dict(
            modfile=os.path.join(
                alphaviz.utils.DATA_PATH,
                'modifications.tsv'
            ),
            aasfile=os.path.join(
                alphaviz.utils.DATA_PATH,
                'amino_acids.tsv'
            ),
            verbose=False,
        )
        self.analysis_software = ""
        self.mz_tol = options.layout[0][0][0]
        self.im_tol = options.layout[0][0][1]
        self.rt_tol = options.layout[0][0][2]
        self.heatmap_x_axis = options.layout[1][0][0]
        self.heatmap_y_axis = options.layout[1][0][1]
        self.heatmap_colormap = options.layout[1][0][2]
        self.heatmap_background_color = options.layout[1][0][3]
        self.heatmap_precursor_size = options.layout[1][0][4]
        self.heatmap_precursor_color = options.layout[1][0][5]
        self.colorscale_qualitative = options.layout[2][0][0]
        self.colorscale_sequential = options.layout[2][0][1]
        self.protein_seq = str()
        self.gene_name = str()
        self.ms1_ms2_frames = dict()
        self.ms1_frame = None
        self.merged_precursor_data = pd.DataFrame()
        self.peptide = dict()
        self.proteins_table = pn.widgets.Tabulator(
            layout='fit_data_table',
            name='Proteins table',
            pagination='remote',
            page_size=5,
            disabled=True,
            height=250,
            show_index=False,
            selectable=1,
            formatters={
                "Protein IDs": {
                    'type': 'link',
                    'urlPrefix':"https://www.uniprot.org/uniprot/",
                    'target':"_blank",
                }
            },
            sizing_mode='stretch_width',
            align='center',
            text_align='center',
            margin=(0, 5, 10, 5)
        )
        self.gene_name_filter = pn.widgets.AutocompleteInput(
            name='Search the protein by its gene name:',
            min_characters=3,
            case_sensitive=False,
            width=350,
            margin=(0, 0, 10, 0),
        )
        self.gene_name_reset = pn.widgets.Button(
            name='\u21bb',
            height=32,
            button_type='default',
            width=50,
            margin=(18, 5, 0, 20),
        )
        self.protein_list_title = pn.pane.Markdown(
            'Load a list of proteins:',
            margin=(10, 5, 0, 20),
        )
        self.protein_list = pn.widgets.FileInput(
            accept='.txt',
            margin=(22, 5, 0, 5),
        )
        self.peptides_table = pn.widgets.Tabulator(
            layout='fit_data_table',
            pagination='remote',
            page_size=8,
            page=1,
            disabled=True,
            height=300,
            show_index=False,
            selectable=1,
            sizing_mode='stretch_width',
            align='center',
            text_align='center',
            margin=(0, 5, 10, 5)
        )
        self.protein_coverage_plot = None
        self.chromatograms_plot = None
        self.heatmap_ms1_plot = None
        self.heatmap_ms2_plot = None
        self.line_plot = None
        self.x_axis_label_mq = pn.widgets.Select(
            name='Select X label:',
            value='Retention time',
            options=['Retention time', 'Ion mobility'],
            width=150,
            align='center'
        )
        self.x_axis_label_diann = pn.widgets.Select(
            name='Select dimension:',
            value='RT dimension',
            options=['RT dimension', 'RT/IM dimension'],
            width=150,
            align='center'
        )
        self.previous_frame = pn.widgets.Button(
            type='default',
            name='\u25c0  Previous frame',
            sizing_mode='stretch_width',
            margin=(10, 30, 30, 30)
        )
        self.next_frame = pn.widgets.Button(
            type='default',
            name='Next frame  \u25b6',
            sizing_mode='stretch_width',
            margin=(10, 30, 30, 30)
        )
        self.plot_overlapped_frames = pn.widgets.Toggle(
            name='Overlap frames',
            button_type='default',
            sizing_mode='stretch_width',
        )
        self.layout = None

    def create_layout(self):
        self.analysis_software = self.data.settings.get('analysis_software')
        self.update_gene_name_filter()
        if self.analysis_software:
            dependances = {
                self.gene_name_reset: [self.reset_protein_table, 'clicks'],
                self.protein_list: [self.filter_protein_table, 'value'],
                self.gene_name_filter: [self.run_after_gene_filter, 'value'],
                self.proteins_table: [self.run_after_protein_selection, 'selection'],
                self.peptides_table: [self.run_after_peptide_selection, 'selection'],
                self.heatmap_x_axis: [self.display_heatmap_spectrum, 'value'],
                self.heatmap_y_axis: [self.display_heatmap_spectrum, 'value'],
                self.heatmap_colormap: [self.display_heatmap_spectrum, 'value'],
                self.heatmap_background_color: [self.display_heatmap_spectrum, 'value'],
                self.heatmap_precursor_size: [self.display_heatmap_spectrum, 'value'],
                self.heatmap_precursor_color: [self.display_heatmap_spectrum, 'value'],
                self.previous_frame: [self.display_previous_frame, 'clicks'],
                self.next_frame: [self.display_next_frame, 'clicks'],
                self.plot_overlapped_frames: [self.display_overlapped_frames, 'value'],
                self.mz_tol: [self.display_line_spectra_plots, 'value'],
                self.im_tol: [self.display_line_spectra_plots, 'value'],
                self.rt_tol: [self.display_line_spectra_plots, 'value'],
                self.x_axis_label_mq: [self.display_line_spectra_plots, 'value'],
                self.x_axis_label_diann: [self.display_elution_profile_plots, 'value'],
                self.colorscale_qualitative: [self.run_after_peptide_selection, 'value'],
                self.colorscale_sequential: [self.run_after_peptide_selection, 'value'],
            }
            for k in dependances.keys():
                k.param.watch(
                    dependances[k][0],
                    dependances[k][1]
                )
        self.dictionary = json.load(open(os.path.join(
            alphaviz.utils.STYLE_PATH,
            'tables_formatting.json',
        )))
        if self.analysis_software == 'maxquant':
            self.proteins_table.value = self.data.mq_protein_groups
            self.proteins_table.formatters = self.dictionary[self.analysis_software]['proteins_table']['formatters']
            self.proteins_table.widths = self.dictionary[self.analysis_software]['proteins_table']['widths']
            self.peptides_table.value = self.data.mq_evidence.iloc[0:0]
            self.peptides_table.widths = self.dictionary[self.analysis_software]['peptides_table']['widths']
            if '(EXP) Seq coverage, %' in self.data.mq_protein_groups.columns:
                self.proteins_table.formatters['(EXP) Seq coverage, %'] = {"type": "progress", "max": 100, "legend": True}

        elif self.analysis_software == 'diann':
            self.proteins_table.value = self.data.diann_proteins
            self.peptides_table.value = self.data.diann_peptides.iloc[0:0]

        # plots
        self.chromatograms_plot = alphaviz.plotting.plot_chrom(
            self.data.raw_data
        )
        if self.analysis_software:
            self.layout = pn.Column(
                pn.Pane(
                    self.chromatograms_plot,
                    config=update_config('Chromatograms'),
                    sizing_mode='stretch_width',
                    margin=(0, 10)
                ),
                pn.Row(
                    self.gene_name_filter,
                    self.gene_name_reset,
                    self.protein_list_title,
                    self.protein_list,
                    margin=(10, 0),
                ),
                pn.panel(
                    f"### Proteins table",
                    align='center',
                    margin=(-10, 10, -5, 10)
                ),
                self.proteins_table,
                pn.panel(
                    f"### Peptides table",
                    align='center',
                    margin=(-10, 10, -5, 10)
                ),
                self.peptides_table,
                self.protein_coverage_plot,
                None, # peptide description
                None, #XIC plot
                pn.Row(
                    None, #Previous frame button
                    None, #Next frame button
                ),
                pn.Row(
                    self.heatmap_ms1_plot,
                    self.heatmap_ms2_plot,
                    sizing_mode='stretch_width',
                    align='center'
                ),
                None, #Overlap frames button
                None, #Summed MS2 spectrum
                margin=(20, 10, 5, 10),
                sizing_mode='stretch_width',
            )
        else:
            self.layout = pn.Column(
                pn.Pane(
                    self.chromatograms_plot,
                    config=update_config('Chromatograms'),
                    sizing_mode='stretch_width',
                    margin=(0, 10)
                ),
                margin=(20, 10, 5, 10),
                sizing_mode='stretch_width',
            )
        return self.layout

    def update_gene_name_filter(self):
        if self.analysis_software == 'maxquant':
            self.gene_name_filter.options = self.data.mq_protein_groups['Gene names'].str.split(';').explode().unique().tolist()
        elif self.analysis_software == 'diann':
            self.gene_name_filter.options = self.data.diann_proteins['Gene names'].str.split(';').explode().unique().tolist()

    def reset_protein_table(self, *args):
        self.proteins_table.loading = True
        self.peptides_table.loading = True
        self.gene_name_filter.value = ''
        if self.analysis_software == 'maxquant':
            self.proteins_table.value = self.data.mq_protein_groups
        elif self.analysis_software == 'diann':
            self.proteins_table.value = self.data.diann_proteins
        self.protein_list.value = b''
        self.proteins_table.selection = []
        self.peptides_table.loading = False
        self.proteins_table.loading = False

    def filter_protein_table(self, *args):
        if self.protein_list.value != b'':
            self.proteins_table.loading = True
            self.peptides_table.loading = True
            self.peptides_table.value = self.data.mq_evidence.iloc[0:0] if self.analysis_software == 'maxquant' else self.data.diann_peptides.iloc[0:0]
            self.proteins_table.selection = []
            predefined_list = []
            for line in StringIO(str(self.protein_list.value, "utf-8")).readlines():
                predefined_list.append(line.strip().upper())
            if predefined_list:
                if self.analysis_software == 'maxquant':
                    filtered_df = alphaviz.preprocessing.filter_df(
                        self.data.mq_protein_groups,
                        pattern='|'.join(predefined_list),
                        column='Gene names',
                        software='maxquant',
                    )
                    if filtered_df.empty:
                        self.proteins_table.value = self.data.mq_protein_groups.iloc[0:0,:]
                    else:
                        self.proteins_table.value = filtered_df
                elif self.analysis_software == 'diann':
                    filtered_df = self.data.diann_proteins[self.data.diann_proteins['Gene names'].isin(predefined_list)]
                    if filtered_df.empty:
                        self.proteins_table.value = self.data.diann_proteins.iloc[0:0,:]
                    else:
                        self.proteins_table.value = filtered_df
            else:
                self.proteins_table.value = self.data.mq_protein_groups if self.analysis_software == 'maxquant' else self.data.diann_proteins
            self.peptides_table.loading = False
            self.proteins_table.loading = False

    def run_after_gene_filter(self, *args):
        self.proteins_table.loading = True
        self.peptides_table.loading = True
        self.proteins_table.selection = []
        if self.analysis_software == 'maxquant':
            self.proteins_table.value = alphaviz.preprocessing.filter_df(
                self.data.mq_protein_groups,
                pattern=self.gene_name_filter.value,
                column='Gene names',
                software='maxquant',
            )
            self.peptides_table.value = self.data.mq_evidence.iloc[0:0]
        elif self.analysis_software == 'diann':
            self.proteins_table.value = alphaviz.preprocessing.filter_df(
                self.data.diann_proteins,
                pattern=self.gene_name_filter.value,
                column='Gene names',
                software='diann',
            )
            self.peptides_table.value = self.data.diann_peptides.iloc[0:0]
        self.peptides_table.loading = False
        self.proteins_table.loading = False

    def run_after_protein_selection(self, *args):
        if self.proteins_table.selection:
            self.peptides_table.loading = True
            self.peptides_table.selection = []
            if self.analysis_software == 'maxquant':
                self.gene_name = self.proteins_table.value.iloc[self.proteins_table.selection[0]]['Gene names']
                self.peptides_table.value = self.data.mq_evidence[self.data.mq_evidence['Gene names'] == self.gene_name]
                if self.peptides_table.value.empty:
                    self.peptides_table.value = self.data.mq_evidence[self.data.mq_evidence['Gene names'].str.contains(self.gene_name)]
                self.curr_protein_ids = [val.replace('CON__', '') if "|" not in val else val.split('|')[1].replace('CON__', '') for val in self.peptides_table.value['Leading razor protein'].sort_values(ascending=False).values][0]
            elif self.analysis_software == 'diann':
                self.gene_name = self.proteins_table.value.iloc[self.proteins_table.selection[0]]['Gene names']
                self.curr_protein_ids = self.proteins_table.value.iloc[self.proteins_table.selection[0]]['Protein IDs']
                self.peptides_table.value = alphaviz.preprocessing.filter_df(
                    self.data.diann_peptides,
                    pattern=self.gene_name,
                    column='Gene names',
                    software='diann',
                )
            self.layout[7:] = [
                None, # peptide description
                None, #XIC plot
                pn.Row(
                    None, #Previous frame button
                    None, #Next frame button
                ),
                pn.Row(
                    None,
                    None,
                    sizing_mode='stretch_width',
                    align='center'
                ),
                None, #Overlap frames button
                None, #Summed MS2 spectrum
            ]
            self.protein_seq = alphaviz.preprocessing.get_aa_seq(
                self.curr_protein_ids,
                self.data.fasta,
            )
            self.protein_coverage_plot = alphaviz.plotting.plot_sequence_coverage(
                self.protein_seq,
                self.gene_name,
                self.peptides_table.value['Modified.Sequence'].tolist() if self.analysis_software == 'diann' else self.peptides_table.value['Modified sequence'].tolist(),
                self.colorscale_qualitative.value,
                self.colorscale_sequential.value,
                r"\[(.*?)\]|\((.*?)\)\)?",
                self.curr_protein_ids
            )
            if not self.protein_coverage_plot and self.analysis_software == 'maxquant':
                curr_protein_ids = sorted(self.peptides_table.value['Proteins'].values[0].split(';'), reverse=True)
                for prot_id in curr_protein_ids:
                    if '|' in prot_id:
                        prot_id = prot_id.split('|')[1]
                    self.protein_seq = alphaviz.preprocessing.get_aa_seq(
                        prot_id,
                        self.data.fasta,
                    )
                    self.protein_coverage_plot = alphaviz.plotting.plot_sequence_coverage(
                        self.protein_seq,
                        self.gene_name,
                        self.peptides_table.value['Modified.Sequence'].tolist() if self.analysis_software == 'diann' else self.peptides_table.value['Modified sequence'].tolist(),
                        self.colorscale_qualitative.value,
                        self.colorscale_sequential.value,
                        r"\[(.*?)\]|\((.*?)\)\)?",
                        prot_id
                    )
                    if self.protein_coverage_plot:
                        self.curr_protein_ids = prot_id
                        break

            self.layout[6] = pn.Pane(
                self.protein_coverage_plot,
                config=update_config(f"{self.gene_name}_coverage_plot"),
                align='center',
                sizing_mode='stretch_width',
            )
            self.peptides_table.loading = False
        else:
            self.peptides_table.loading = True
            self.peptides_table.selection = []
            self.peptides_table.value = self.data.mq_evidence.iloc[0:0] if self.analysis_software == 'maxquant' else self.data.diann_peptides.iloc[0:0]
            self.layout[6] = None
            self.layout[7:] = [
                None, # peptide description
                None, #XIC plot
                pn.Row(
                    None, #Previous frame button
                    None, #Next frame button
                ),
                pn.Row(
                    None,
                    None,
                    sizing_mode='stretch_width',
                    align='center'
                ),
                None, #Overlap frames button
                None, #Summed MS2 spectrum
            ]
            self.peptides_table.loading = False

    def run_after_peptide_selection(self, *args):
        if self.proteins_table.selection:
            self.peptides_table.loading = True
            if self.peptides_table.selection:
                self.protein_coverage_plot = alphaviz.plotting.plot_sequence_coverage(
                    self.protein_seq,
                    self.gene_name,
                    [self.peptides_table.value.iloc[self.peptides_table.selection[0]]['Modified.Sequence']] if self.analysis_software == 'diann' else [self.peptides_table.value.iloc[self.peptides_table.selection[0]]['Modified sequence']],
                    self.colorscale_qualitative.value,
                    self.colorscale_sequential.value,
                    r"\[(.*?)\]|\((.*?)\)\)?",
                    self.curr_protein_ids
                )
                self.layout[6] = pn.Pane(
                    self.protein_coverage_plot,
                    config=update_config(f"{self.gene_name}_coverage_plot"),
                    align='center',
                    sizing_mode='stretch_width',
                )
                self.scan_number = [int(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['MS/MS scan number'])]
                if 'dda' in self.data.raw_data.acquisition_mode:
                    pasef_ids = [int(pasef_id) for pasef_id in self.data.mq_all_peptides[self.data.mq_all_peptides['MS/MS scan number'].isin(self.scan_number)]['Pasef MS/MS IDs'].values[0]]
                    precursors = self.data.raw_data.fragment_frames[self.data.raw_data.fragment_frames.index.isin(pasef_ids)].copy()
                    # quick fix the AlphaTims's bug with the differences in the Frames in raw_data.fragment_frames table for .d and .hdf files
                    if self.data.ms_file_name.value.split('.')[-1] == 'hdf':
                        precursors.loc[:, 'Frame'] -= 1
                    self.merged_precursor_data = pd.merge(
                        precursors, self.data.raw_data.precursors[self.data.raw_data.precursors.Id.isin(precursors.Precursor.values)],
                        left_on='Precursor',
                        right_on='Id'
                    )
                    self.merged_precursor_data['Frame_Prec'] = list(zip(self.merged_precursor_data.Frame, self.merged_precursor_data.Precursor))
                    self.ms1_ms2_frames = dict(zip(self.merged_precursor_data.Parent, self.merged_precursor_data.Frame_Prec))
                    self.current_frame = list(self.ms1_ms2_frames.keys())[0]
                    self.display_line_spectra_plots()
                    self.display_heatmap_spectrum()
                else:
                    self.ms2_frame = self.data.raw_data.fragment_frames[self.data.raw_data.fragment_frames.index.isin(self.scan_number)].Frame.values[0]
                    if self.data.ms_file_name.value.split('.')[-1] == 'hdf':
                        self.ms2_frame -= 1
                    self.ms1_frame = self.data.raw_data.frames[(self.data.raw_data.frames.MsMsType == 0) & (self.data.raw_data.frames.Id < self.ms2_frame)].iloc[-1, 0]
                    self.peptide = {
                        "sequence":
                        self.peptides_table.value.iloc[self.peptides_table.selection[0]]['Sequence_AP_mod'],
                        "charge":
                        self.peptides_table.value.iloc[self.peptides_table.selection[0]]['Charge'],
                        "im":
                        self.peptides_table.value.iloc[self.peptides_table.selection[0]]['IM'],
                        "rt": self.peptides_table.value.iloc[self.peptides_table.selection[0]]['RT'] * 60
                    }
                    self.peptide['mz'] = alphaviz.utils.calculate_mz(
                        prec_mass=alphaviz.utils.get_precmass(
                            alphaviz.utils.parse(self.peptide['sequence']),
                            self.mass_dict
                        ),
                        charge=self.peptide['charge']
                    )
                    self.display_elution_profile_plots()
                    self.display_heatmap_spectrum()
            else:
                self.peptides_table.selection = []
                self.layout[7:] = [
                    None, # peptide description
                    None, #XIC plot
                    pn.Row(
                        None, #Previous frame button
                        None, #Next frame button
                    ),
                    pn.Row(
                        None,
                        None,
                        sizing_mode='stretch_width',
                        align='center'
                    ),
                    None, #Overlap frames button
                    None, #Summed MS2 spectrum
                ]
            self.peptides_table.loading = False

    def display_line_spectra_plots(self, *args):
        if self.analysis_software == 'maxquant' and not self.merged_precursor_data.empty:
            try:
                self.layout[8][1].loading = True
            except IndexError:
                pass
            mz_tol_value = self.mz_tol.value
            prec_mono_mz = self.merged_precursor_data.MonoisotopicMz.median()
            prec_mono_low_mz = prec_mono_mz / (1 + mz_tol_value / 10**6)
            prec_mono_high_mz = prec_mono_mz * (1 + mz_tol_value / 10**6)
            if self.x_axis_label_mq.value == 'Retention time':
                one_over_k0 = float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['1/K0'])
                one_over_k0_low, one_over_k0_high = one_over_k0 - self.im_tol.value, one_over_k0 + self.im_tol.value
                precursor_indices = self.data.raw_data[
                    :,
                    one_over_k0_low : one_over_k0_high,
                    :,
                    prec_mono_low_mz : prec_mono_high_mz,
                    'raw'
                ]
            else:
                precursor_indices = self.data.raw_data[
                    :,
                    :,
                    :,
                    prec_mono_low_mz : prec_mono_high_mz,
                    'raw'
                ]
            self.layout[7] = pn.panel(
                f"## The selected peptide: m/z: {round(float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['m/z']), 3)}, charge: {float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['Charge'])}, 1/K0: {round(float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['1/K0']), 3)}, andromeda score: {round(float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['Andromeda score']), 1)}.",
                css_classes=['main-part'],
                sizing_mode='stretch_width',
                align='center',
                margin=(0, 10, 0, -10)
            )
            conversion_dict = {
                'Retention time': 'rt',
                'Ion mobility': 'mobility'
            }
            self.layout[8] = pn.Row(
                self.x_axis_label_mq,
                pn.Pane(
                    alphaviz.plotting.plot_line(
                        self.data.raw_data,
                        precursor_indices,
                        conversion_dict[self.x_axis_label_mq.value],
                    ),
                    sizing_mode='stretch_width',
                    config=update_config('Extracted Ion Chromatogram'),
                    loading=False
                ),
                sizing_mode='stretch_width',
                margin=(5, 10, 0, 10)
            )
        else:
            self.display_elution_profile_plots()

    def display_elution_profile_plots(self, *args):
        if self.analysis_software == 'diann' and self.peptide:
            self.layout[7] = pn.panel(
                f"## The selected peptide: m/z: {round(self.peptide['mz'], 3)}, charge: {self.peptide['charge']}, 1/K0: {round(self.peptide['im'], 3)}, Quantity.Quality score: {round(float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['Quantity.Quality']), 2)}.",
                css_classes=['main-part'],
                sizing_mode='stretch_width',
                align='center',
                margin=(0, 10, 0, -10)
            )

            if self.x_axis_label_diann.value == 'RT dimension':
                self.layout[8] = pn.Row(
                    self.x_axis_label_diann,
                    pn.Pane(
                        alphaviz.plotting.plot_elution_profile(
                            self.data.raw_data,
                            self.peptide,
                            self.mass_dict,
                            mz_tol=self.mz_tol.value,
                            rt_tol=self.rt_tol.value,
                            im_tol=self.im_tol.value,
                            title=f"Precursor and fragment elution profiles of {self.peptides_table.value.iloc[self.peptides_table.selection[0]]['Modified.Sequence']} in RT dimension ({self.peptide['rt'] / 60:.2f} min)",
                            colorscale_qualitative=self.colorscale_qualitative.value,
                            colorscale_sequential=self.colorscale_sequential.value,
                        ),
                        sizing_mode='stretch_width',
                        config=update_config('Precursor&fragment elution profile plot in RT dimension'),
                        loading=False,
                    ),
                    margin=(5, 10, 0, 10)
                )
            else:
                self.layout[8] = pn.Row(
                    self.x_axis_label_diann,
                    pn.pane.HoloViews(
                        alphaviz.plotting.plot_elution_profile_heatmap(
                            self.data.raw_data,
                            self.peptide,
                            self.mass_dict,
                            mz_tol=self.mz_tol.value,
                            rt_tol=self.rt_tol.value,
                            im_tol=self.im_tol.value,
                            n_cols=8,
                            width=180,
                            height=180,
                            colormap=self.heatmap_colormap.value,
                            background_color=self.heatmap_background_color.value,
                        ),
                        sizing_mode='stretch_width',
                        linked_axes=True,
                        loading=False,
                    ),
                    sizing_mode='stretch_width',
                    margin=(5, 10, 0, 10)
                )

    def display_heatmap_spectrum(self, *args):
        if self.ms1_ms2_frames or self.ms1_frame:
            if self.analysis_software == 'maxquant':
                ms1_frame = self.current_frame
                ms2_frame = self.ms1_ms2_frames[self.current_frame][0]
                mz = float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['m/z'])
                im = float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['1/K0'])
            elif self.analysis_software == 'diann':
                ms1_frame = self.ms1_frame
                ms2_frame = self.ms2_frame
                mz = self.peptide['mz']
                im = self.peptide['im']
            data_ms1 = self.data.raw_data[ms1_frame].copy()
            try:
                self.heatmap_ms1_plot = alphaviz.plotting.plot_heatmap(
                    data_ms1,
                    mz=mz,
                    im=im,
                    x_axis_label=self.heatmap_x_axis.value,
                    y_axis_label=self.heatmap_y_axis.value,
                    title=f'MS1 frame(s) #{ms1_frame}',
                    colormap=self.heatmap_colormap.value,
                    background_color=self.heatmap_background_color.value,
                    precursor_size=self.heatmap_precursor_size.value,
                    precursor_color=self.heatmap_precursor_color.value,
                    width=570,
                    height=450,
                    margin=(0, 10, 10, 0),
                    # shared_axes=True
                )
                data_ms2 = self.data.raw_data[ms2_frame].copy()
                self.heatmap_ms2_plot = alphaviz.plotting.plot_heatmap(
                    data_ms2,
                    x_axis_label=self.heatmap_x_axis.value,
                    y_axis_label=self.heatmap_y_axis.value,
                    title=f'MS2 frame(s) #{ms2_frame}',
                    colormap=self.heatmap_colormap.value,
                    background_color=self.heatmap_background_color.value,
                    width=570,
                    height=450,
                    margin=(0, 10, 10, 0),
                    # shared_axes=True
                )

                self.layout[10][0] = pn.pane.HoloViews(
                    self.heatmap_ms1_plot,
                    margin=(15, 0, 0, 0),
                    linked_axes=False,
                    loading=False
                )
                self.layout[10][1] = pn.pane.HoloViews(
                    self.heatmap_ms2_plot,
                    margin=(15, 0, 0, 0),
                    linked_axes=False,
                    loading=False
                )
            except ValueError:
                print('The x- and y-axis of the heatmaps should be different.')
            if self.analysis_software == 'diann':
                if self.x_axis_label_diann.value == 'RT/IM dimension':
                    self.display_elution_profile_plots()
            if self.analysis_software == 'maxquant':
                data_ions = alphaviz.preprocessing.get_mq_ms2_scan_data(
                    self.data.mq_msms,
                    self.scan_number[0],
                    self.data.raw_data,
                    self.ms1_ms2_frames[self.current_frame][1]
                )
                self.ms_spectra_plot = alphaviz.plotting.plot_mass_spectra(
                    data_ions,
                    title=f'MS2 spectrum for Precursor: {self.ms1_ms2_frames[self.current_frame][1]}',
                    sequence=self.peptides_table.value.iloc[self.peptides_table.selection[0]]['Sequence']
                )
                self.layout[9][0] = self.previous_frame
                self.layout[9][1] = self.next_frame
                self.layout[11] = self.plot_overlapped_frames
                self.layout[12] = pn.Pane(
                    self.ms_spectra_plot,
                    config=update_config('Combined MS2 spectrum'),
                    margin=(10, 0, 0, 0),
                    sizing_mode='stretch_width',
                    loading=False,
                    height=600
                )

    def display_previous_frame(self, *args):
        try:
            self.layout[10][0].loading = True
            self.layout[10][1].loading = True
            self.layout[12].loading = True
        except IndexError:
            pass
        self.plot_overlapped_frames.value = False
        current_frame_index = list(self.ms1_ms2_frames.keys()).index(self.current_frame)
        if current_frame_index == 0:
            self.current_frame = list(self.ms1_ms2_frames.keys())[-1]
        else:
            self.current_frame = list(self.ms1_ms2_frames.keys())[current_frame_index - 1]
        self.display_heatmap_spectrum()

    def display_next_frame(self, *args):
        try:
            self.layout[10][0].loading = True
            self.layout[10][1].loading = True
            self.layout[12].loading = True
        except IndexError:
            pass
        self.plot_overlapped_frames.value = False
        current_frame_index = list(self.ms1_ms2_frames.keys()).index(self.current_frame)
        if current_frame_index == len(self.ms1_ms2_frames.keys())-1:
            self.current_frame = list(self.ms1_ms2_frames.keys())[0]
        else:
            self.current_frame = list(self.ms1_ms2_frames.keys())[current_frame_index + 1]
        self.display_heatmap_spectrum()

    def display_overlapped_frames(self, *args):
        try:
            self.layout[10][0].loading = True
            self.layout[10][1].loading = True
            self.layout[12].loading = True
        except IndexError:
            pass
        if self.plot_overlapped_frames.value == True:
            mz = float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['m/z'])
            im = float(self.peptides_table.value.iloc[self.peptides_table.selection[0]]['1/K0'])
            try:
                self.heatmap_ms1_plot = alphaviz.plotting.plot_heatmap(
                    self.data.raw_data[list(self.ms1_ms2_frames.keys())],
                    mz=mz,
                    im=im,
                    x_axis_label=self.heatmap_x_axis.value,
                    y_axis_label=self.heatmap_y_axis.value,
                    title=f'MS1 frame(s) #{list(self.ms1_ms2_frames.keys())}',
                    colormap=self.heatmap_colormap.value,
                    background_color=self.heatmap_background_color.value,
                    width=570,
                    height=450,
                    margin=(0, 10, 10, 0),
                )
                self.heatmap_ms2_plot = alphaviz.plotting.plot_heatmap(
                    self.data.raw_data[[val[0] for val in self.ms1_ms2_frames.values()]],
                    x_axis_label=self.heatmap_x_axis.value,
                    y_axis_label=self.heatmap_y_axis.value,
                    title=f'MS2 frame(s) #{[val[0] for val in self.ms1_ms2_frames.values()]}',
                    colormap=self.heatmap_colormap.value,
                    background_color=self.heatmap_background_color.value,
                    width=570,
                    height=450,
                    margin=(0, 10, 10, 0),
                )
                self.layout[10][0] = pn.pane.HoloViews(
                    self.heatmap_ms1_plot,
                    margin=(15, 0, 0, 0),
                    linked_axes=False,
                    loading=False
                )
                self.layout[10][1] = pn.pane.HoloViews(
                    self.heatmap_ms2_plot,
                    margin=(15, 0, 0, 0),
                    linked_axes=False,
                    loading=False
                )
            except ValueError:
                print('The x- and y-axis of the heatmaps should be different.')
            self.layout[12] = None
        else:
            self.display_heatmap_spectrum()

class QCTab(object):

    def __init__(self, data, options):
        self.name = "Quality Control"
        self.data = data
        self.mz_tol = options.layout[0][0][0]
        self.layout_qc = None
        self.analysis_software = self.data.settings.get('analysis_software')

    def create_layout(self):
        dependances = {
            self.mz_tol: [self.update_tol_range_mass_density, 'value'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1]
            )

        experiment = self.data.ms_file_name.value.split('.')[0]
        if self.analysis_software == 'maxquant':
            uncalb_mass_dens_plot = alphaviz.plotting.plot_mass_error(
                self.data.mq_evidence,
                'm/z',
                'Uncalibrated mass error [ppm]',
                'Uncalibrated mass density plot',
                self.mz_tol.value,
            )
            calb_mass_dens_plot = alphaviz.plotting.plot_mass_error(
                self.data.mq_evidence,
                'm/z',
                'Mass error [ppm]',
                'Calibrated mass density plot',
            )
            peptide_per_protein_distr = alphaviz.plotting.plot_pept_per_protein_barplot(
                self.data.mq_protein_groups,
                '(EXP) # peptides',
                'Peptides per protein',
            )
            peptide_mz_distr = alphaviz.plotting.plot_peptide_distr(
                self.data.mq_evidence,
                'm/z',
                'Peptide m/z distribution',
            )
            peptide_length_distr = alphaviz.plotting.plot_peptide_distr(
                self.data.mq_evidence,
                'Length',
                'Peptide length distribution',
            )

            self.layout_qc = pn.Column(
                pn.widgets.Tabulator(
                    self.data.mq_summary,
                    sizing_mode='stretch_width',
                    layout='fit_data_table',
                    name='Overview table',
                    selection=list(self.data.mq_summary[self.data.mq_summary['Raw file'].str.contains(experiment)].index),
                    row_height=40,
                    disabled=True,
                    height=200,
                    show_index=False,
                ),
                pn.panel(
                    f"## Quality control of the entire sample",
                    align='center',
                    margin=(15, 10, -5, 10)
                ),
                pn.Row(
                    pn.Pane(
                        uncalb_mass_dens_plot,
                        config=update_config('Uncalibrated mass density plot'),
                    ),
                    pn.Pane(
                        calb_mass_dens_plot,
                        config=update_config('Calibrated mass density plot'),
                    ),
                    align='center'
                ),
                pn.Row(
                    peptide_per_protein_distr,
                    peptide_mz_distr,
                    peptide_length_distr,
                    align='center',
                ),
                margin=(0, 10, 5, 10),
                sizing_mode='stretch_width',
                align='start',
            )
        elif self.analysis_software == 'diann':
            peptide_per_protein_distr = alphaviz.plotting.plot_pept_per_protein_barplot(
                self.data.diann_proteins,
                '(EXP) # peptides',
                'Peptides per protein',
            )
            peptide_charge_distr = alphaviz.plotting.plot_peptide_distr(
                self.data.diann_peptides,
                'Charge',
                'Peptide charge distribution'
            )
            peptide_length_distr = alphaviz.plotting.plot_peptide_distr(
                self.data.diann_peptides,
                'Length',
                'Peptide length distribution'
            )
            self.layout_qc = pn.Column(
                pn.widgets.Tabulator(
                    self.data.diann_statist,
                    sizing_mode='stretch_width',
                    layout='fit_data_table',
                    name='Overview table',
                    selection=list(self.data.diann_statist[self.data.diann_statist['File.Name'].str.contains(experiment)].index),
                    row_height=40,
                    disabled=True,
                    show_index=False,
                ),
                pn.panel(
                    f"## Quality control of the entire sample",
                    align='center',
                    margin=(15, 10, -5, 10)
                ),
                pn.Row(
                    peptide_per_protein_distr,
                    peptide_charge_distr,
                    peptide_length_distr,
                    align='center',
                ),
                margin=(0, 10, 5, 10),
                sizing_mode='stretch_width',
                align='start',
            )
        else:
            self.layout_qc = pn.pane.Markdown(
                'To use this functionality, load the output data from any supported software analysis tool.',
                margin=(5, 0, 0, 10),
            )
        return self.layout_qc

    def update_tol_range_mass_density(self, *args):
        if self.analysis_software == 'maxquant' and self.layout_qc:
            self.layout_qc[2][0] = pn.Pane(
                alphaviz.plotting.plot_mass_error(
                    self.data.mq_evidence,
                    'm/z',
                    'Uncalibrated mass error [ppm]',
                    'Uncalibrated mass density plot',
                    self.mz_tol.value,
                ),
                config=update_config('Uncalibrated mass density plot'),
            )


class TargetModeTab(object):
    def __init__(self, data, options):
        self.name = "Targeted Mode"
        self.data = data
        self.mz_tol = options.layout[0][0][0]
        self.im_tol = options.layout[0][0][1]
        self.rt_tol = options.layout[0][0][2]
        self.heatmap_x_axis = options.layout[1][0][0]
        self.heatmap_y_axis = options.layout[1][0][1]
        self.heatmap_colormap = options.layout[1][0][2]
        self.heatmap_background_color = options.layout[1][0][3]
        self.heatmap_precursor_size = options.layout[1][0][4]
        self.heatmap_precursor_color = options.layout[1][0][5]
        self.colorscale_qualitative = options.layout[2][0][0]
        self.colorscale_sequential = options.layout[2][0][1]
        self.layout_target_mode = None
        self.analysis_software = self.data.settings.get('analysis_software')
        self.targeted_peptides_table = pn.widgets.Tabulator(
            value=pd.DataFrame(
                columns=['name','sequence', 'charge', 'im', 'rt']
            ),
            widths={'index': 70},
            layout='fit_columns',
            selectable=1,
            height=250,
            show_index=True,
            width=570,
            margin=(25, 12, 10, 18)
        )
        self.peptides_count = pn.widgets.IntInput(
            name='Add N empty row(s)',
            value=0,
            step=1,
            start=0,
            end=1000
        )
        self.peptides_table_text = pn.pane.Markdown(
            'Load a table of targeted peptides:',
            margin=(5, 0, 0, 10),
        )
        self.peptides_table_file = pn.widgets.FileInput(
            accept='.tsv,.csv,.txt',
            margin=(-10,0,0,10)
        )
        self.mass_dict = alphaviz.utils.get_mass_dict(
            modfile=os.path.join(
                alphaviz.utils.DATA_PATH,
                'modifications.tsv'
            ),
            aasfile=os.path.join(
                alphaviz.utils.DATA_PATH,
                'amino_acids.tsv'
            ),
            verbose=False,
        )
        self.clear_peptides_table_button = pn.widgets.Button(
            name='Clear table',
            button_type='default',
            width=60,
            margin=(25, 12, 10, 18)
        )

    def create_layout(self):
        experiment = self.data.ms_file_name.value.split('.')[0]

        dependances = {
            self.peptides_table_file: [self.read_peptides_table, 'value'],
            self.targeted_peptides_table: [self.visualize_elution_plots, ['selection', 'value']],
            self.peptides_count: [self.update_row_count, 'value'],
            self.heatmap_colormap: [self.visualize_elution_plots, 'value'],
            self.heatmap_background_color: [self.visualize_elution_plots, 'value'],
            self.mz_tol: [self.visualize_elution_plots, 'value'],
            self.im_tol: [self.visualize_elution_plots, 'value'],
            self.rt_tol: [self.visualize_elution_plots, 'value'],
            self.colorscale_qualitative: [self.visualize_elution_plots, 'value'],
            self.colorscale_sequential: [self.visualize_elution_plots, 'value'],
            self.clear_peptides_table_button: [
            self.clear_peptide_table, 'clicks'
            ]
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1]
            )

        if 'dia' in self.data.raw_data.acquisition_mode:
            self.layout_target_mode = pn.Column(
                pn.Row(
                    pn.Column(
                        self.peptides_count,
                        self.peptides_table_text,
                        self.peptides_table_file,
                    ),
                    self.targeted_peptides_table,
                    self.clear_peptides_table_button,
                ),
                None,
                None,
                margin=(15, 10, 5, 10),
                sizing_mode='stretch_width',
                align='start',
            )
        else:
            self.layout_target_mode = pn.Column(
                pn.pane.Markdown(
                    'To use this functionality please load DIA data.',
                    margin=(5, 0, 0, 10),
                ),
                None,
                None,
            )
        return self.layout_target_mode

    def clear_peptide_table(self, *args):
        if not self.targeted_peptides_table.value.empty:
            self.targeted_peptides_table.value = pd.DataFrame(
                columns=['name', 'sequence', 'charge', 'im', 'rt'],
            )

    def update_row_count(self, *args):
        if self.targeted_peptides_table.value.empty:
            self.targeted_peptides_table.value = pd.DataFrame(
                columns=['name', 'sequence', 'charge', 'im', 'rt'],
                index=range(self.peptides_count.value),
            )
        else:
            self.targeted_peptides_table.value = self.targeted_peptides_table.value.append(
                pd.DataFrame(
                    columns=self.targeted_peptides_table.value.columns,
                    index=range(self.peptides_count.value),
                ),
                ignore_index=True
            )


    def read_peptides_table(self, *args):
        file_ext = os.path.splitext(self.peptides_table_file.filename)[-1]
        if file_ext=='.csv':
            sep=';'
        else:
            sep='\t'
        self.targeted_peptides_table.selection = []
        self.targeted_peptides_table.value = pd.read_csv(
            StringIO(str(self.peptides_table_file.value, "utf-8")),
            sep=sep
        )

    def visualize_elution_plots(self, *args):
        if 'dia' in self.data.raw_data.acquisition_mode:
            if self.targeted_peptides_table.selection:
                try:
                    peptide = self.targeted_peptides_table.value.iloc[self.targeted_peptides_table.selection[0]].to_dict()
                except IndexError:
                    peptide = {}
                if peptide and not any(pd.isna(val) for val in peptide.values()):
                    self.targeted_peptides_table.loading = True
                    try:
                        peptide['charge'] = int(peptide['charge'])
                        for val in ['im', 'rt']:
                            peptide[val] = float(peptide[val])
                        for val in ['name', 'sequence']:
                            peptide[val] = str(peptide[val])
                        peptide['mz'] = alphaviz.utils.calculate_mz(
                            prec_mass=alphaviz.utils.get_precmass(
                                alphaviz.utils.parse(peptide['sequence']),
                                self.mass_dict
                            ),
                            charge=peptide['charge']
                        )
                    except:
                        print('The current peptide cannot be loaded.')
                    else:
                        peptide['rt'] *= 60 # to convert to seconds
                        try:
                            self.layout_target_mode[1] = pn.Pane(
                                alphaviz.plotting.plot_elution_profile(
                                    self.data.raw_data,
                                    peptide,
                                    self.mass_dict,
                                    mz_tol=self.mz_tol.value,
                                    rt_tol=self.rt_tol.value,
                                    im_tol=self.im_tol.value,
                                    title=f"Precursor fragment elution profiles of {peptide['name']}({peptide['sequence']}) in RT and RT/IM dimensions ({peptide['rt'] / 60: .2f} min)",
                                    colorscale_qualitative=self.colorscale_qualitative.value,
                                    colorscale_sequential=self.colorscale_sequential.value,
                                    height=500,
                                ),
                                sizing_mode='stretch_width',
                                config=update_config('Precursor/fragments elution profile plot'),
                                loading=False,
                            )
                        except:
                            self.layout_target_mode[1] = None
                        try:
                            self.layout_target_mode[2] = pn.pane.HoloViews(
                                alphaviz.plotting.plot_elution_profile_heatmap(
                                    self.data.raw_data,
                                    peptide,
                                    self.mass_dict,
                                    mz_tol=self.mz_tol.value,
                                    rt_tol=self.rt_tol.value,
                                    im_tol=self.im_tol.value,
                                    n_cols=8,
                                    width=180,
                                    height=180,
                                    colormap=self.heatmap_colormap.value,
                                    background_color=self.heatmap_background_color.value,
                                ),
                                sizing_mode='stretch_width',
                                linked_axes=True,
                                loading=False,
                                align='center',
                            )
                        except AttibuteError:
                            self.layout_target_mode[2] = None
                    finally:
                        self.targeted_peptides_table.loading = False
                else:
                    self.layout_target_mode[1], self.layout_target_mode[2] = None, None
            else:
                self.layout_target_mode[1], self.layout_target_mode[2] = None, None


class GUI(object):
    # TODO: move to alphabase and docstring

    def __init__(
        self,
        name,
        github_url,
        run_in_background=False,
        automatic_close=True,
    ):
        self.name = name
        self.tab_counter = 0
        self.header = HeaderWidget(
            name,
            alphaviz.utils.IMG_PATH,
            github_url
        )
        self.layout = pn.Column(
            self.header.create_layout(),
            sizing_mode='stretch_width',
            min_width=1270
        )
        self.run_in_background = run_in_background
        self.automatic_close = automatic_close

    def start_server(self, run_in_background=False):
        if self.automatic_close:
            bokeh_ws_handler = bokeh.server.views.ws.WSHandler
            self.bokeh_server_open = bokeh_ws_handler.open
            bokeh_ws_handler.open = self.__open_browser_tab(
                self.bokeh_server_open
            )
            self.bokeh_server_on_close = bokeh_ws_handler.on_close
            bokeh_ws_handler.on_close = self.__close_browser_tab(
                self.bokeh_server_on_close
            )
        self.server = self.layout.show(threaded=True, title=self.name)
        if not run_in_background:
            self.server.join()
        elif not self.run_in_background:
            self.server.join()

    def __open_browser_tab(self, func):
        def wrapper(*args, **kwargs):
            self.tab_counter += 1
            return func(*args, **kwargs)
        return wrapper

    def __close_browser_tab(self, func):
        def wrapper(*args, **kwargs):
            self.tab_counter -= 1
            return_value = func(*args, **kwargs)
            if self.tab_counter == 0:
                self.stop_server()
            return return_value
        return wrapper

    def stop_server(self):
        logging.info("Stopping server...")
        self.server.stop()
        if self.automatic_close:
            bokeh_ws_handler = bokeh.server.views.ws.WSHandler
            bokeh_ws_handler.open = self.bokeh_server_open
            bokeh_ws_handler.on_close = self.bokeh_server_on_close


class AlphaVizGUI(GUI):
    # TODO: docstring

    def __init__(self, start_server=False):
        super().__init__(
            name=f"AlphaViz {alphaviz.__version__}",
            github_url='https://github.com/MannLabs/alphaviz',
        )
        # TODO: include main widget in base gui?
        self.project_description = """### AlphaViz is an automated visualization pipeline to link the identifications found in the analysis with the original raw MS data and easily assess their individual quality or the overall quality of whole samples."""
        self.manual_path = os.path.join(
            alphaviz.utils.DOCS_PATH,
            "alphaviz_tutorial.pdf"
        )
        self.main_widget = MainWidget(
            self.project_description,
            self.manual_path
        )

        # ERROR/WARNING MESSAGES
        self.error_message_upload = "The selected file can't be uploaded. Please check the instructions for data uploading."

        self.data = DataImportWidget()
        self.options = OptionsWidget(self.data)
        self.options.add_option(ToleranceOptionsWidget().create_layout())
        self.options.add_option(HeatmapOptionsWidget().create_layout())
        self.options.add_option(ColorscaleOptionsWidget().create_layout())
        self.tabs = TabsWidget(self.data, self.options)
        self.layout += [
            self.main_widget.create_layout(),
            self.data.create_layout(),
            self.options.get_layout(),
            self.tabs.create_layout(
                [
                    ('Main View', pn.panel("Blank")),
                    ('Quality Control', pn.panel("Blank")),
                    ('Targeted Mode', pn.panel("Blank"))
                ]
            ),
        ]
        if start_server:
            self.start_server()


def run():
    init_panel()
    AlphaVizGUI(start_server=True)


if __name__ == '__main__':
    run()
