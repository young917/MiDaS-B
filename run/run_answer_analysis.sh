dataset=("email-Enron-full" "email-Eu-full" "contact-high-school" "contact-primary-school" "NDC-classes-full" "NDC-substances-full" "tags-ask-ubuntu" "tags-math-sx" "threads-ask-ubuntu" "coauth-MAG-Geology-full" "coauth-MAG-History-full")
spset=("0.1" "0.3" "0.5" "0.7" "0.9")

cd ../
for data in ${dataset[@]}
do
    # Evaluation
    ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir answer --algo_opt "intersection,densification,sizewcc,clusteringcoef" --accuracy 100 
    cd src
    python calculation_helper.py --dataset ${data} --algorithm answer --effdiam --overlapness  --accuracy 100
    cd ..
    for sp in ${spset[@]}
    do
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir answer --samplingportion ${sp} 
        
        # For threads-ask-ubuntu, coauth-MAG-Geology-full, and coauth-MAG-History-full, 
        # run `src/preprocess_sv.py --algoname answer --dataname ${data}`
        # and then run `run_sv_answer.m` before the next line

        cd src
        python calculation_helper.py --dataset ${data} --algorithm answer --svdist --samplingportion ${sp} 
        cd ..
    done
done