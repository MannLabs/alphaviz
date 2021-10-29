#!bash

# Initial cleanup
rm -rf dist
rm -rf build
FILE=alphaviz.pkg
if test -f "$FILE"; then
  rm alphaviz.pkg
fi
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n alphavizinstaller python=3.8 -y
conda activate alphavizinstaller

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_macos_gui
pip install "../../dist/alphaviz-0.0.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.2
pyinstaller ../pyinstaller/alphaviz.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../alphaviz/data/*.fasta dist/alphaviz/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/alphaviz/Contents/Resources
cp ../logos/alpha_logo.icns dist/alphaviz/Contents/Resources
mv dist/alphaviz_gui dist/alphaviz/Contents/MacOS
cp Info.plist dist/alphaviz/Contents
cp alphaviz_terminal dist/alphaviz/Contents/MacOS
cp ../../LICENSE.txt Resources/LICENSE.txt
cp ../logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/alphaviz --identifier de.mpg.biochem.alphaviz.app --version 0.0.1 --install-location /Applications/alphaviz.app --scripts scripts alphaviz.pkg
productbuild --distribution distribution.xml --resources Resources --package-path alphaviz.pkg dist/alphaviz_gui_installer_macos.pkg
