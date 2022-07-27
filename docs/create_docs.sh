rm -rf build
conda env remove -n alphavizdocs
conda create -n alphavizdocs python=3.8 -y
# conda create -n alphavizinstaller python=3.8
conda activate alphavizdocs
# call conda install git -y
# call pip install 'git+https://github.com/MannLabs/alphaviz.git#egg=alphaviz[gui]' --use-feature=2020-resolver
# brew install freetype
pip install '../.[stable,gui-stable]'
make html
conda deactivate
