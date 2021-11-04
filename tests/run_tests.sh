conda activate alphaviz
# python -m unittest test_cli
# python -m unittest test_gui
python -m pytest test_io.py
python -m pytest test_preprocessing.py
conda deactivate
