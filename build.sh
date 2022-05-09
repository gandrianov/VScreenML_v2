$PYTHON setup.py install

chmod +x bin/calculate_features.py
chmod +x bin/predict_class.py
chmod +x bin/pyrosetta_minimizer.py
chmod +x bin/mol2genparams.py

cp bin/calculate_features.py $PREFIX/bin/
cp bin/predict_class.py $PREFIX/bin/
cp bin/pyrosetta_minimizer.py $PREFIX/bin/
cp bin/mol2genparams.py $PREFIX/bin/


