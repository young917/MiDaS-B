dataset=("email-Enron-full" "email-Eu-full" "contact-high-school" "contact-primary-school" "NDC-classes-full" "NDC-substances-full" "tags-ask-ubuntu" "tags-math-sx" "threads-ask-ubuntu" "coauth-MAG-Geology-full" "coauth-MAG-History-full")
spset=("0.1" "0.3" "0.5" "0.7" "0.9")
alphaset=("0.00" "0.25" "0.50" "1.00" "2.00")
betaset=("-1.00" "-0.50" "-0.25" "0.00" "0.25" "0.50" "1.00")

cd ../
for repeat_index in 1 2 3
do
    for data in ${dataset[@]}
    do
        for alpha in ${alphaset[@]}
        do
            for beta in ${betaset[@]}
            do
                # Sampling
                ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm essz --algo_opt add_global_deg_min --alpha ${alpha} --beta ${beta} --repeat ${repeat_index}
                # Score
                ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir essz/add_global_deg_min_${alpha}_${beta} --algo_opt "intersection,densification,sizewcc" --accuracy 500 --repeat ${repeat_index}
                cd src
                python calculation_avg_helper.py --dataset $data --algorithm essz/add_global_deg_min_${alpha}_${beta} --size --repeat ${repeat_index} --accuracy 500
                python score_function.py --dataset ${data} --datasetdir ../dataset/ --algorithm essz/add_global_deg_min_${alpha}_${beta} --repeat ${repeat_index} --accuracy 500
                cd ..  
                # Evaluation
                ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir essz/add_global_deg_min_${alpha}_${beta} --algo_opt "clusteringcoef" --accuracy 100 --repeat ${repeat_index}
                cd src
                python calculation_helper.py --dataset ${data} --algorithm essz/add_global_deg_min_${alpha}_${beta} --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
                cd ..
                for sp in ${spset[@]}
                do
                    ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir essz/add_global_deg_min_${alpha}_${beta}  --samplingportion ${sp} --repeat ${repeat_index}
                    # For threads-ask-ubuntu, coauth-MAG-Geology-full, and coauth-MAG-History-full, 
                    # run  `src/preprocess_sv.py --dataname ${data} --algoname essz/add_global_deg_min_${alpha}_${beta} --repeat_str ${repeat_index}` 
                    # and then run `run_sv_midasB_oracle.m` before the next line
                    cd src
                    python calculation_helper.py --dataset ${data} --algorithm essz/add_global_deg_min_${alpha}_${beta} --svdist --samplingportion ${sp} --repeat ${repeat_index}
                    cd ..
                done
            done
        done
    done
done

cd analyze/
python evaluation_table.py --search_name midasB_oracle --repeat "1,2,3"
# run `analyze/GridSearch.ipynb`
cd results/
xargs -a midasB_oracle/search_result.txt cp -r -t midasB_oracle/
