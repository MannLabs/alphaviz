conda create -n alphaviz_pip_test python=3.8 -y
conda activate alphaviz_pip_test
pip install "alphaviz[stable,stable-gui]"
alphaviz
conda deactivate
