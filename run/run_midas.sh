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
            # Score
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir es/add_global_deg_min_${alpha} --algo_opt "intersection,densification,sizewcc" --accuracy 500 --repeat ${repeat_index} --recalculate
            cd src
            python calculation_avg_helper.py --dataset $data --algorithm es/add_global_deg_min_${alpha} --size --repeat ${repeat_index} --accuracy 500 --recalculate
            python score_function.py --dataset ${data} --datasetdir ../dataset/ --algorithm es/add_global_deg_min_${alpha} --repeat ${repeat_index} --accuracy 500
            cd ..  
        done
    done
done

# Find best one based on score
cd ./analyze/
python BestScoreSearch_midas.py
cd ..

cd ./results/
xargs -a midas/search_result.txt cp -r -t midas/

# # Evaluation
for repeat_index in 1 2 3
do
    for data in ${dataset[@]}
    do
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir midas --algo_opt "clusteringcoef" --accuracy 100 --repeat ${repeat_index}
        cd src
        python calculation_helper.py --dataset ${data} --algorithm midas --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
        cd ..
        for sp in ${spset[@]}
        do
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir midas  --samplingportion ${sp} --repeat ${repeat_index}
            
            # For threads-ask-ubuntu, coauth-MAG-Geology-full, and coauth-MAG-History-full, 
            # run  `src/preprocess_sv.py --dataname ${data} --algoname midas --repeat_str ${repeat_index}` 
            # and then run `run_sv_midas.m` before the next line
            
            cd src
            python calculation_helper.py --dataset ${data} --algorithm midas --svdist --samplingportion ${sp} --repeat ${repeat_index}
            cd ..
        done
    done
done