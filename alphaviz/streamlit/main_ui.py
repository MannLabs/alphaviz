import streamlit as st

from PIL import Image
import os

import alphaviz

from alphaviz.streamlit import (
    raw_viz, tims_viz
)

_this_file = __file__
_this_directory = os.path.dirname(_this_file)
LOGO_PATH = os.path.join(_this_directory, 'logos', 'alpha_logo.png')
ICON_PATH = os.path.join(_this_directory, 'logos', 'alpha_logo.ico')
image = Image.open(LOGO_PATH)
icon = Image.open(ICON_PATH)

st.set_page_config(
    page_title=f"AlphaViz {alphaviz.__version__}",
    # page_icon=icon,
    layout="wide",
)

hide_streamlit_style = """
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>

"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.sidebar.image(image, width = 300)
st.sidebar.code(f"AlphaViz {alphaviz.__version__}")

sidebar = {
    'Raw': raw_viz.show,
    'TimsTOF': tims_viz.show,
}

menu = st.sidebar.radio("", list(sidebar.keys()))

if menu:
    sidebar[menu]()

link = f'[AlphaViz on GitHub]({alphaviz.__github__})'
st.sidebar.markdown(link, unsafe_allow_html=True)