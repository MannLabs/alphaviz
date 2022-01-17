import os
import logging
import pandas as pd
from io import StringIO

import alphatims.bruker
import alphatims.utils

# visualization
import panel as pn
import bokeh.server.views.ws
from bokeh.models.widgets.tables import NumberFormatter
import plotly.express as px
import holoviews as hv

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

    def create_layout(self):
        self.layout = pn.Row(
            self.project_description,
            pn.layout.HSpacer(width=500),
            self.manual,
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
            name='Specify a folder with Bruker raw files:',
            placeholder='Enter the full path to the folder with Bruker .d folders/.hdf files',
            width=900,
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.ms_file_name = pn.widgets.Select(
            name='Select a raw file',
            size=10,
            width=900,
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.path_output_folder = pn.widgets.TextInput(
            name='Specify an analysis folder:',
            placeholder='Enter the full path to the analysis folder',
            width=900,
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.path_fasta_file = pn.widgets.TextInput(
            # TODO: remove the fixed fasta file before release
            value='/Users/eugeniavoytik/copied/Bruker/MaxQuant_output_tables/20210413_TIMS03_EVO03_PaSk_MA_HeLa_200ng_S1-A1_1_24848.d/txt/human.fasta',
            name='Specify a fasta file:',
            placeholder='Enter the full path to the used for the analysis fasta file',
            width=900,
            sizing_mode='stretch_width',
            margin=(15, 15, 15, 15)
        )
        self.is_prediction = pn.widgets.Checkbox(
            name='Activate the prediction',
            margin=(5, 0, 5, 15),
        )
        # UPLOAD DATA
        self.upload_button = pn.widgets.Button(
            name='Upload Data',
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
            sizing_mode='stretch_width',
            # object='test warning message',
            margin=(30, 15, 5, 15),
        )

    def create_layout(self):
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
            # height=470,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        self.path_raw_folder.param.watch(
            self.update_file_names,
            'value'
        )
        self.ms_file_name.param.watch(
            self.update_output_folder_and_fasta,
            'value'
        )
        self.upload_button.param.watch(
            self.load_data,
            'clicks'
        )

        return self.layout

    def update_file_names(self, *args):
        self.ms_file_name.options = sorted(
            alphaviz.io.get_filenames_from_directory(
                self.path_raw_folder.value,
                ['d', 'hdf']
            )
        )

    def update_output_folder_and_fasta(self, *args):
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

    def load_data(self, *args):
        alphatims.utils.set_progress_callback(self.upload_progress)
        self.upload_progress.value = 0
        self.raw_data = alphatims.bruker.TimsTOF(
            os.path.join(
                self.path_raw_folder.value,
                self.ms_file_name.value
            ),
        )
        alphatims.utils.set_progress_callback(True)
        # TODO: to change it just by changing the active=True parameter for the
        self.upload_progress = pn.indicators.Progress(
            active=True,
            bar_color='light',
            width=250,
            align='center',
            margin=(-10, 0, 30, 0)
        )
        self.layout[0][2][1] = self.upload_progress

        # read a fasta file
        self.fasta = alphaviz.io.read_fasta(
            self.path_fasta_file.value
        )

        # read analysis output files (MQ, DIA-NN, etc.)
        ## check all files in the analysis output folder
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
                raise FileNotFoundError('The MQ output files necessary for the visualization are not found.')
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
                raise FileNotFoundError('The DIA-NN output files necessary for the visualization are not found.')

        self.trigger_dependancy()
        self.upload_progress.active = False


class OptionsWidget(object):

    def __init__(self, data):
        self.info = pn.panel('Will be updated soon.')
        self.data = data
        self.layout = pn.Card(
            title='Options',
            collapsed=True,
            sizing_mode='stretch_width',
            # height=225,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

    def get_layout(self, *args):
        return self.layout

    def add_option(self, option):
        self.layout.append(option)


class HeatmapOptionsWidget(object):

    def __init__(self):
        # self.plot1_title = pn.pane.Markdown(
        #     '#### Axis for Heatmap',
        #     margin=(10, 0, -5, 0),
        #     align='center'
        # )
        self.plot1_x_axis = pn.widgets.Select(
            name='X axis',
            value='m/z, Th',
            options=['m/z, Th', 'Inversed IM, V·s·cm\u207B\u00B2', 'RT, min'],
            width=180,
            margin=(20, 20, 20, 20),
        )
        self.plot1_y_axis = pn.widgets.Select(
            name='Y axis',
            value='Inversed IM, V·s·cm\u207B\u00B2',
            options=['m/z, Th', 'Inversed IM, V·s·cm\u207B\u00B2', 'RT, min'],
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
            options=['black', 'white', 'red', 'yellow', 'green', 'blue'],
            width=180,
            margin=(20, 20, 20, 10),
        )
        self.precursor_target_size = pn.widgets.IntInput(
            name='Precursor target size',
            value=15,
            step=5,
            start=0,
            end=100,
            width=180,
            margin=(20, 20, 20, 10),
        )
        self.precursor_target_color = pn.widgets.Select(
            name='Precursor target color',
            value='blue',
            options=['black', 'white', 'red', 'yellow', 'green', 'blue'],
            width=180,
            margin=(20, 20, 20, 10),
        )

    def create_layout(self, *args):
        layout = pn.Card(
            # self.plot1_title,
            pn.Row(
                self.plot1_x_axis,
                self.plot1_y_axis,
                self.heatmap_color,
                self.heatmap_background,
                self.precursor_target_size,
                self.precursor_target_color
            ),
            title='Heatmap Options',
            collapsed=False,
            sizing_mode='stretch_width',
            # height=225,
            margin=(15, 8, 0, 8),
            css_classes=['background']
        )
        return layout


class XicOptionsWidget(object):

    def __init__(self):
        self.mz_tolerance = pn.widgets.FloatInput(
            name='M/z Tolerance (ppm)',
            value=10,
            step=5,
            start=0,
            end=1000,
            width=100
        )
        self.mz_tolerance_units = pn.widgets.Select(
            name='M/z Tolerance Units',
            value='ppm',
            options=['ppm', 'Da'],
            width=130
        )
        self.im_tolerance = pn.widgets.FloatInput(
            name='IM Tolerance (ppm)',
            value=0.05,
            step=0.1,
            start=0,
            end=2,
            width=100
        )
        self.rt_tolerance = pn.widgets.FloatInput(
            name='RT Tolerance (ppm)',
            value=30,
            step=5,
            start=0,
            end=1000,
            width=100
        )

    def create_layout(self, *args):
        layout = pn.Card(
            pn.Row(
                self.mz_tolerance,
                self.mz_tolerance_units,
                self.im_tolerance,
                self.rt_tolerance
            ),
            title='XIC Options',
            collapsed=False,
            sizing_mode='stretch_width',
            # height=225,
            margin=(15, 8, 15, 8),
            css_classes=['background']
        )
        return layout

class CoveragePlotOptionsWidget(object):

    def __init__(self):
        self.peptides_qualitative_colorscale = pn.widgets.Select(
            name='Peptides qualitative color scale',
            value='Alphabet',
            options=list(set([each['y'][0] for each in px.colors.qualitative.swatches()['data']])),
            width=180,
            margin=(20, 20, 20, 10),
        )
        self.peptides_sequential_colorscale = pn.widgets.Select(
            name='Peptides sequential color scale',
            value='Viridis',
            options=list(set([each['y'][0] for each in px.colors.sequential.swatches()['data']])),
            width=180,
            margin=(20, 20, 20, 10),
        )

    def create_layout(self, *args):
        layout = pn.Card(
            pn.Row(
                self.peptides_qualitative_colorscale,
                self.peptides_sequential_colorscale,
            ),
            title='Coverage Plot Options',
            collapsed=False,
            sizing_mode='stretch_width',
            # height=225,
            margin=(15, 8, 15, 8),
            css_classes=['background']
        )
        return layout


class TabsWidget(object):

    def __init__(self, data, options):
        self.info = pn.panel('Will be updated soon.')
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
        if (self.data.mq_evidence is not None or self.data.diann_proteins is not None) and self.data.raw_data is not None:
            self.layout = pn.Tabs(
                tabs_location='above',
                margin=(30, 10, 5, 8),
                sizing_mode='stretch_width'
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
        self.name = "Main View"
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
        self.analysis_software = self.data.settings.get('analysis_software')
        self.heatmap_x_axis = options.layout[0][0][0]
        self.heatmap_y_axis = options.layout[0][0][1]
        self.heatmap_colormap = options.layout[0][0][2]
        self.heatmap_background_color = options.layout[0][0][3]
        self.heatmap_precursor_size = options.layout[0][0][4]
        self.heatmap_precursor_color = options.layout[0][0][5]
        self.mz_tol = options.layout[1][0][0]
        self.mz_tol_units = options.layout[1][0][1]
        self.im_tol = options.layout[1][0][2]
        self.rt_tol = options.layout[1][0][3]
        self.colorscale_qualitative = options.layout[2][0][0]
        self.colorscale_sequential = options.layout[2][0][1]
        self.protein_seq = str()
        self.gene_name = str()
        self.ms1_ms2_frames = dict()
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
            name='Search a protein by a gene name:',
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
        self.x_axis_label = pn.widgets.Select(
            name='Select X label',
            value='rt',
            options=['rt', 'mobility'],
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
        self.update_gene_name_filter()

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
            self.mz_tol_units: [self.display_line_spectra_plots, 'value'],
            self.im_tol: [self.display_line_spectra_plots, 'value'],
            self.rt_tol: [self.display_line_spectra_plots, 'value'],
            self.x_axis_label: [self.display_line_spectra_plots, 'value'],
            self.colorscale_qualitative: [self.run_after_protein_selection, 'value'],
            self.colorscale_sequential: [self.run_after_protein_selection, 'value'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1]
            )

        self.dictionary = {
            'maxquant': {
                'peptides_table': {
                    'formatters': {
                        'Acetylation (N-term)': {
                            'type': 'tickCross',
                            'allowTruthy': True
                        },
                        'Oxidation (M)': {
                            'type': 'tickCross',
                            'allowTruthy': True
                        },
        #                 'Mass': NumberFormatter(format='0,0.000'),
        #                 'm/z': NumberFormatter(format='0,0.000'),
        #                 '1/K0': NumberFormatter(format='0,0.000'),
        #                 'Intensity': NumberFormatter(format='0,0'),
        #                 'MS/MS scan number': NumberFormatter(format='0,0'),
        #                 'Andromeda score':  NumberFormatter(format='0,0.0'),
                    },
                    'widths': {
                        'Sequence': 220,
                        'Proteins': 200,
                        'MS/MS scan number': 100,
                        'Oxidation (M)': 130,
                    },
                },
                'proteins_table': {
                    'formatters': {
                        '(EXP) Seq coverage, %': {
                            'type': 'progress',
                            'max': 100,
                            'legend': True
                        },
                        'Protein names': {
                            'type': "textarea"
                        },
                    },
                    'widths': {
                        'Protein IDs': 230,
                        'Protein names': 350,
                        'Sequence lengths': 150,
                    },
                }
            },
            'diann': {
                'peptides_table': {},
                'proteins_table': {}
            }
        }
        if self.analysis_software == 'maxquant':
            self.proteins_table.value = self.data.mq_protein_groups
            self.proteins_table.formatters = self.dictionary[self.analysis_software]['proteins_table']['formatters']
            self.proteins_table.widths = self.dictionary[self.analysis_software]['proteins_table']['widths']
            self.peptides_table.value = self.data.mq_evidence.iloc[0:0]
            self.peptides_table.formatters = self.dictionary[self.analysis_software]['peptides_table']['formatters']
            self.peptides_table.widths = self.dictionary[self.analysis_software]['peptides_table']['widths']
        elif self.analysis_software == 'diann':
            self.proteins_table.value = self.data.diann_proteins
            self.peptides_table.value = self.data.diann_peptides.iloc[0:0]

        # plots
        self.chromatograms_plot = alphaviz.plotting.plot_chrom(
            self.data.raw_data
        )

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
        # print('inside filter_protein_table')
        if self.protein_list.value != b'':
            # print('inside if self.protein_list.value != b''')
            self.proteins_table.loading = True
            self.peptides_table.loading = True
            self.peptides_table.value = self.data.mq_evidence.iloc[0:0] if self.analysis_software == 'maxquant' else self.data.diann_peptides.iloc[0:0]
            self.proteins_table.selection = []
            # print(self.proteins_table.selection)
            predefined_list = []
            for line in StringIO(str(self.protein_list.value, "utf-8")).readlines():
                predefined_list.append(line.strip().upper())
            if self.analysis_software == 'maxquant':
                # print(self.proteins_table.value.shape)
                self.proteins_table.value = alphaviz.preprocessing.filter_df(
                    self.data.mq_protein_groups,
                    pattern='|'.join(predefined_list),
                    column='Gene names',
                    software='maxquant',
                )
                # print(self.proteins_table.value.shape)
            elif self.analysis_software == 'diann':
                self.proteins_table.value = self.data.diann_proteins[self.data.diann_proteins['Gene names'].isin(predefined_list)]
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
                curr_protein_ids = self.proteins_table.value.iloc[self.proteins_table.selection[0]]['Protein IDs']
                # print(self.gene_name, curr_protein_ids)
                self.peptides_table.value = alphaviz.preprocessing.filter_df(
                    self.data.mq_evidence,
                    pattern=self.gene_name.replace(';', '|'), # use the regex | character to try to match each of the substrings in the genes separated by ; in the "Gene names" column
                    column='Gene names',
                    software='maxquant',
                )
                # print(self.peptides_table.value)
            elif self.analysis_software == 'diann':
                # print('11')
                self.gene_name = self.proteins_table.value.iloc[self.proteins_table.selection[0]]['Gene naames']
                # print('22')
                curr_protein_ids = self.proteins_table.value.iloc[self.proteins_table.selection[0]]['Protein IDs']
                # print(self.gene_name, curr_protein_ids)
                self.peptides_table.value = alphaviz.preprocessing.filter_df(
                    self.data.diann_peptides,
                    pattern=self.gene_name,
                    column='Gene names',
                    software='diann',
                )
                # print('3')
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
            # print('4')
            self.protein_seq = alphaviz.preprocessing.get_aa_seq(
                curr_protein_ids,
                self.data.fasta,
            )
            # print('5')
            self.protein_coverage_plot = alphaviz.plotting.plot_sequence_coverage(
                self.protein_seq,
                self.gene_name,
                self.peptides_table.value['Modified.Sequence'].tolist() if self.analysis_software == 'diann' else self.peptides_table.value['Modified sequence'].tolist(),
                self.colorscale_qualitative.value,
                self.colorscale_sequential.value,
                r"\[([^]]+)\]|\((\w+)\)"
            )
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
            self.protein_coverage_plot = alphaviz.plotting.plot_sequence_coverage(
                self.protein_seq,
                self.gene_name,
                self.peptides_table.selected_dataframe['Modified.Sequence'].tolist() if self.analysis_software == 'diann' else self.peptides_table.selected_dataframe['Modified sequence'].tolist(),
                self.colorscale_qualitative.value,
                self.colorscale_sequential.value,
                r"\[([^]]+)\]|\((\w+)\)"
            )
            self.layout[6] = pn.Pane(
                self.protein_coverage_plot,
                config=update_config(f"{self.gene_name}_coverage_plot"),
                align='center',
                sizing_mode='stretch_width',
            )

            if self.peptides_table.selected_dataframe.shape[0] == 1:
                self.scan_number = [int(scan) for scan in self.peptides_table.selected_dataframe['MS/MS scan number'].tolist()]
                if 'dda' in self.data.raw_data.acquisition_mode:
                    pasef_ids = [int(pasef_id) for pasef_id in self.data.mq_all_peptides[self.data.mq_all_peptides['MS/MS scan number'].isin(self.scan_number)]['Pasef MS/MS IDs'].values[0]]
                    precursors = self.data.raw_data.fragment_frames[self.data.raw_data.fragment_frames.index.isin(pasef_ids)]
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
                    self.ms1_frame = self.data.raw_data.frames[(self.data.raw_data.frames.MsMsType == 0) & (self.data.raw_data.frames.Id < self.ms2_frame)].iloc[-1, 0]
                    self.peptide = {
                        "sequence": self.peptides_table.selected_dataframe['Sequence_AP_mod'].values[0],
                        "charge":
                        self.peptides_table.selected_dataframe['Charge'].values[0],
                        "im": self.peptides_table.selected_dataframe['IM'].values[0],
                        "rt": self.peptides_table.selected_dataframe['RT'].values[0] * 60
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
        if self.analysis_software == 'maxquant':
            try:
                self.layout[8][1].loading = True
            except IndexError:
                pass
            mz_tol_value = self.mz_tol.value
            prec_mono_mz = self.merged_precursor_data.MonoisotopicMz.median()
            if self.mz_tol_units.value == 'ppm':
                prec_mono_low_mz = prec_mono_mz / (1 + mz_tol_value / 10**6)
                prec_mono_high_mz = prec_mono_mz * (1 + mz_tol_value / 10**6)
            else:
                prec_mono_low_mz = prec_mono_mz - mz_tol_value
                prec_mono_high_mz = prec_mono_mz + mz_tol_value
            if self.x_axis_label.value == 'rt':
                one_over_k0 = float(self.peptides_table.selected_dataframe['1/K0'].values[0])
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
                f"## The selected peptide: m/z: {round(float(self.peptides_table.selected_dataframe['m/z'].values[0]), 3)}, charge: {float(self.peptides_table.selected_dataframe['Charge'].values[0])}, 1/K0: {round(float(self.peptides_table.selected_dataframe['1/K0'].values[0]), 3)}, andromeda score: {round(float(self.peptides_table.selected_dataframe['Andromeda score'].values[0]), 1)}.",
                css_classes=['main-part'],
                sizing_mode='stretch_width',
                align='center',
                margin=(0, 10, 0, -10)
            )

            self.layout[8] = pn.Row(
                self.x_axis_label,
                pn.Pane(
                    alphaviz.plotting.plot_line(
                        self.data.raw_data,
                        precursor_indices,
                        self.x_axis_label.value,
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
        self.layout[7] = pn.panel(
            f"## The selected peptide: m/z: {round(self.peptide['mz'], 3)}, charge: {self.peptide['charge']}, 1/K0: {round(self.peptide['im'], 3)}, Quantity.Quality score: {round(float(self.peptides_table.selected_dataframe['Quantity.Quality'].values[0]), 2)}.",
            css_classes=['main-part'],
            sizing_mode='stretch_width',
            align='center',
            margin=(0, 10, 0, -10)
        )

        if self.x_axis_label.value == 'rt':
            self.layout[8] = pn.Row(
                self.x_axis_label,
                pn.Pane(
                    alphaviz.plotting.plot_elution_profile(
                        self.data.raw_data,
                        self.peptide,
                        self.mass_dict,
                        mz_tol=self.mz_tol.value,
                        rt_tol=self.rt_tol.value,
                        im_tol=self.im_tol.value,
                        title=f"Precursor/fragments elution profile of {self.peptides_table.selected_dataframe['Modified.Sequence'].values[0]} in RT dimension ({self.peptide['rt'] / 60: .2f} min)"
                    ),
                    sizing_mode='stretch_width',
                    config=update_config('Precursor/fragments elution profile plot'),
                    loading=False,
                ),
                # sizing_mode='stretch_width',
                margin=(5, 10, 0, 10)
            )
        else:
            self.layout[8] = pn.Row(
                self.x_axis_label,
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
                        height=180
                        # title=f"Precursor/fragments elution profile of {self.peptides_table.selected_dataframe['Modified.Sequence'].values[0]} in RT dimension ({self.peptide['rt'] / 60: .2f} min)"
                    ),
                    sizing_mode='stretch_width',
                    linked_axes=True,
                    config=update_config('Precursor/fragments elution profile plot'),
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
                mz = float(self.peptides_table.selected_dataframe['m/z'].values[0])
                im = float(self.peptides_table.selected_dataframe['1/K0'].values[0])
            elif self.analysis_software == 'diann':
                ms1_frame = self.ms1_frame
                ms2_frame = self.ms2_frame
                mz = self.peptide['mz']
                im = self.peptide['im']

            data_ms1 = self.data.raw_data[ms1_frame].copy()
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
                margin=(0, 10, 10, 0)
                # shared_axes=True
            )
            data_ms2 = self.data.raw_data[ms2_frame].copy()
            self.heatmap_ms2_plot = alphaviz.plotting.plot_heatmap(
                data_ms2,
                mz=mz,
                im=im,
                x_axis_label=self.heatmap_x_axis.value,
                y_axis_label=self.heatmap_y_axis.value,
                title=f'MS2 frame(s) #{ms2_frame}',
                colormap=self.heatmap_colormap.value,
                background_color=self.heatmap_background_color.value,
                precursor_size=self.heatmap_precursor_size.value,
                precursor_color=self.heatmap_precursor_color.value,
                width=570,
                height=450,
                # shared_axes=True
            )

            self.layout[10][0] = pn.pane.HoloViews(
                self.heatmap_ms1_plot,
                margin=(15, 0, 0, 0),
                # sizing_mode='stretch_width',
                linked_axes=False,
                loading=False
            )
            self.layout[10][1] = pn.pane.HoloViews(
                self.heatmap_ms2_plot,
                margin=(15, 0, 0, 0),
                # sizing_mode='stretch_width',
                linked_axes=False,
                loading=False
            )

            if self.analysis_software == 'maxquant':
                data_ions = alphaviz.preprocessing.get_mq_ms2_scan_data(
                    self.data.mq_msms,
                    self.scan_number[0],
                    self.data.raw_data,
                    self.ms1_ms2_frames[self.current_frame][1]
                )
                # print(self.layout)
                self.ms_spectra_plot = alphaviz.plotting.plot_mass_spectra(
                    data_ions,
                    title=f'MS2 spectrum for Precursor: {self.ms1_ms2_frames[self.current_frame][1]}',
                    sequence=self.peptides_table.selected_dataframe['Sequence'].values[0]
                )
                self.layout[9][0] = self.previous_frame
                self.layout[9][1] = self.next_frame
                self.layout[11] = self.plot_overlapped_frames
                self.layout[12] = pn.Pane(
                    self.ms_spectra_plot,
                    config=update_config('Combined MS2 spectrum'),
                    margin=(10, 0, 0, 0),
                    sizing_mode='stretch_width',
                    loading=False
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
            self.heatmap_ms1_plot = alphaviz.plotting.plot_heatmap(
                self.data.raw_data[list(self.ms1_ms2_frames.keys())],
                mz=float(self.peptides_table.selected_dataframe['m/z'].values[0]),
                im=float(self.peptides_table.selected_dataframe['1/K0'].values[0]),
                x_axis_label=self.heatmap_x_axis.value,
                y_axis_label=self.heatmap_y_axis.value,
                title=f'MS1 frame(s) #{list(self.ms1_ms2_frames.keys())}',
                colormap=self.heatmap_colormap.value,
                background_color=self.heatmap_background_color.value,
                width=570,
                height=450,
            )
            self.heatmap_ms2_plot = alphaviz.plotting.plot_heatmap(
                self.data.raw_data[[val[0] for val in self.ms1_ms2_frames.values()]],
                mz=float(self.peptides_table.selected_dataframe['m/z'].values[0]),
                im=float(self.peptides_table.selected_dataframe['1/K0'].values[0]),
                x_axis_label=self.heatmap_x_axis.value,
                y_axis_label=self.heatmap_y_axis.value,
                title=f'MS2 frame(s) #{[val[0] for val in self.ms1_ms2_frames.values()]}',
                colormap=self.heatmap_colormap.value,
                background_color=self.heatmap_background_color.value,
                width=570,
                height=450,
            )
            self.layout[10][0] = pn.pane.HoloViews(
                self.heatmap_ms1_plot,
                margin=(15, 0, 0, 0),
                # sizing_mode='stretch_width',
                linked_axes=False,
                loading=False
            )
            self.layout[10][1] = pn.pane.HoloViews(
                self.heatmap_ms2_plot,
                margin=(15, 0, 0, 0),
                # sizing_mode='stretch_width',
                linked_axes=False,
                loading=False
            )
        else:
            self.display_heatmap_spectrum()

class QCTab(object):

    def __init__(self, data, options):
        self.name = "Quality Control"
        self.data = data
        self.layout_qc = None
        self.analysis_software = self.data.settings.get('analysis_software')

    def create_layout(self):
        experiment = self.data.ms_file_name.value.split('.')[0]
        # print(self.data.mq_protein_groups.columns)
        # print(self.data.mq_protein_groups.head())
        if self.analysis_software == 'maxquant':
            uncalb_mass_dens_plot = alphaviz.plotting.plot_mass_error(
                self.data.mq_evidence,
                'm/z',
                'Uncalibrated mass error [ppm]',
                'Uncalibrated mass density plot'
            )
            calb_mass_dens_plot = alphaviz.plotting.plot_mass_error(
                self.data.mq_evidence,
                'm/z',
                'Mass error [ppm]',
                'Calibrated mass density plot'
            )
            peptide_per_protein_distr = alphaviz.plotting.plot_pept_per_protein_barplot(
                self.data.mq_protein_groups,
                '(EXP) # peptides',
                'Peptides per protein',
            )
            peptide_mz_distr = alphaviz.plotting.plot_peptide_distr(
                self.data.mq_evidence,
                'm/z',
                'Peptide m/z distribution'
            )
            peptide_length_distr = alphaviz.plotting.plot_peptide_distr(
                self.data.mq_evidence,
                'Length',
                'Peptide length distribution'
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
                    # margin=(0, 0, 0, 50)
                ),
                margin=(0, 10, 5, 10),
                sizing_mode='stretch_width',
                align='start',
            )
        else:
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
            # experiment = self.data.ms_file_name.value.split('.')[0]
            self.layout_qc = pn.Column(
                pn.widgets.Tabulator(
                    self.data.diann_statist,
                    sizing_mode='stretch_width',
                    layout='fit_data_table',
                    name='Overview table',
                    selection=list(self.data.diann_statist[self.data.diann_statist['File.Name'].str.contains(experiment)].index),
                    row_height=40,
                    disabled=True,
                    # height=100,
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
                    # margin=(0, 0, 0, 50)
                ),
                margin=(0, 10, 5, 10),
                sizing_mode='stretch_width',
                align='start',
            )
        return self.layout_qc


class TargetModeTab(object):

    def __init__(self, data, options):
        self.name = "Targeted Mode"
        self.data = data
        self.layout_target_mode = None
        self.analysis_software = self.data.settings.get('analysis_software')

    def create_layout(self):
        experiment = self.data.ms_file_name.value.split('.')[0]
        if 'dia' in self.data.raw_data.acquisition_mode:
            # print("INSIDE")
            self.layout_target_mode = pn.Column(
                # pn.widgets.Tabulator(
                #     self.data.mq_summary,
                #     sizing_mode='stretch_width',
                #     layout='fit_data_table',
                #     name='Overview table',
                #     selection=list(self.data.mq_summary[self.data.mq_summary['Raw file'].str.contains(experiment)].index),
                #     row_height=40,
                #     disabled=True,
                #     height=200,
                #     show_index=False,
                # ),
                # pn.panel(
                #     f"## Quality control of the entire sample",
                #     align='center',
                #     margin=(15, 10, -5, 10)
                # ),
                # pn.Row(
                #     pn.Pane(
                #         uncalb_mass_dens_plot,
                #         config=update_config('Uncalibrated mass density plot'),
                #     ),
                #     pn.Pane(
                #         calb_mass_dens_plot,
                #         config=update_config('Calibrated mass density plot'),
                #     ),
                #     align='center'
                # ),
                # pn.Row(
                #     peptide_per_protein_distr,
                #     peptide_mz_distr,
                #     peptide_length_distr,
                #     align='center',
                #     # margin=(0, 0, 0, 50)
                # ),
                margin=(0, 10, 5, 10),
                sizing_mode='stretch_width',
                align='start',
            )

        return self.layout_target_mode

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
        self.options.add_option(HeatmapOptionsWidget().create_layout())
        self.options.add_option(XicOptionsWidget().create_layout())
        self.options.add_option(CoveragePlotOptionsWidget().create_layout())
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
