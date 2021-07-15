#cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO ..

#CC=/usr/bin/clang CXX=/usr/bin/clang++ cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO ..



python3 analyses/perform_analyses.py --threads 16 --methods EXPERIMENTAL  --samples 100 --sample-base ../data/SCER --executable ./build/tsxCount --output runs_weismann_new --max-time 60m

#python3 analyses/perform_analyses.py --threads 8 --methods EXPERIMENTAL  --samples 100 --sample-base ./data/SCER --executable ./build/tsxCount --output runs_test --max-time 60m

/usr/bin/time --verbose ./build/tsxCount --input=./data/SCER.100.fastq --mode=EXPERIMENTAL --threads=8 --l=27 --check --checkabort
/usr/bin/time --verbose ./build/tsxCount --input=../data/SCER.100.fastq --mode=EXPERIMENTAL --threads=8 --l=27 --check --checkabort