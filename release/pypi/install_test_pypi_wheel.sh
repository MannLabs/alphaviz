conda create -n alphaviz_pip_test python=3.8 -y
conda activate alphaviz_pip_test
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple "alphaviz[stable,gui-stable]"
alphaviz
conda deactivate
