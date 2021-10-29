conda create -n alphaviz python=3.8 -y
conda activate alphaviz
pip install -e '../.[stable,development-stable]'
alphaviz
conda deactivate
