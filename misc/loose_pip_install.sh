conda create -n alphaviz python=3.8 -y
conda activate alphaviz
pip install -e '../.[development,gui]'
alphaviz
conda deactivate
