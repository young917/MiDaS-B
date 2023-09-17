dataset=("email-Enron-full" "email-Eu-full" "contact-high-school" "contact-primary-school" "NDC-classes-full" "NDC-substances-full" "tags-ask-ubuntu" "tags-math-sx" "threads-ask-ubuntu" "coauth-MAG-Geology-full" "coauth-MAG-History-full")
spset=("0.1" "0.3" "0.5" "0.7" "0.9")
alphaset=("0.0000" "0.2500" "0.5000" "1.0000" "2.0000")

cd ../
for repeat_index in 1 2 3
do
    for data in ${dataset[@]}
    do
        for alpha in ${alphaset[@]}
        do
            # Sampling
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm es --algo_opt add_global_deg_min --alpha ${alpha} --repeat ${repeat_index}
            # Evaluation
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir es/add_global_deg_min_${alpha} --algo_opt "intersection,densification,sizewcc,clusteringcoef" --accuracy 100 --repeat ${repeat_index}
            cd src
            python calculation_helper.py --dataset ${data} --algorithm es/add_global_deg_min_${alpha} --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
            cd ..
            for sp in ${spset[@]}
            do
                ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir es/add_global_deg_min_${alpha}  --samplingportion ${sp} --repeat ${repeat_index}
                # For threads-ask-ubuntu, coauth-MAG-Geology-full, and coauth-MAG-History-full, 
                # run  `src/preprocess_sv.py --dataname ${data} --algoname es/add_global_deg_min_${alpha} --repeat_str ${repeat_index}` 
                # and then run `run_sv_midas_oracle.m` before the next line
                cd src
                python calculation_helper.py --dataset ${data} --algorithm es/add_global_deg_min_${alpha} --svdist --samplingportion ${sp} --repeat ${repeat_index}
                cd ..
            done
        done
    done
done

cd analyze/
python evaluation_table.py --search_name midas_oracle --repeat "1,2,3"
python GridSearch.py --search_name midas_oracle
cd ../results/
xargs -a midas_oracle/search_result.txt cp -r -t midas_oracle/

