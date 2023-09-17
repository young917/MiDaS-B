# dataset=("email-Enron-full" "email-Eu-full" "contact-high-school" "contact-primary-school" "NDC-classes-full" "NDC-substances-full" "tags-ask-ubuntu" "tags-math-sx" "threads-ask-ubuntu" "coauth-MAG-Geology-full" "coauth-MAG-History-full")
spset=("0.1" "0.3" "0.5" "0.7" "0.9")

dataset=("email-Eu-full" "contact-primary-school" "NDC-substances-full" "tags-ask-ubuntu") # "coauth-MAG-History-full")

cd ../
# for repeat_index in 1 2 3
for repeat_index in 1
do
    for data in ${dataset[@]}
    do
        # Sampling
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm ns --algo_opt add_global_deg --alpha 0.0000 --repeat ${repeat_index} # RNS
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm ns --algo_opt add_global_deg --alpha 1.0000  --repeat ${repeat_index} # RDN
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm rw --algo_opt rw_c --maxlength 1.0 --restart 0.3 --repeat ${repeat_index} # RW
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm ff --algo_opt ff_c --p 0.51 --q 0.20 --repeat ${repeat_index} # FF
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm ns --algo_opt add_ordered --alpha -1.0 --repeat ${repeat_index} # ONS
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm es --algo_opt add_global_deg_min --alpha 0.0000 --repeat ${repeat_index} # RHS
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm tihs --repeat ${repeat_index} # TIHS
        # Evaluation
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir ns/add_global_deg_0.0000 --algo_opt "intersection,densification,clusteringcoef,sizewcc" --accuracy 100 --repeat ${repeat_index}
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir ns/add_global_deg_1.0000 --algo_opt "intersection,densification,clusteringcoef,sizewcc" --accuracy 100 --repeat ${repeat_index}
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir rw/rw_c_1.00_0.30 --algo_opt "intersection,densification,clusteringcoef,sizewcc" --accuracy 100 --repeat ${repeat_index}
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir ff/ff_c_0.51_0.20 --algo_opt "intersection,densification,clusteringcoef,sizewcc" --accuracy 100 --repeat ${repeat_index}
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir ns/add_ordered --algo_opt "intersection,densification,clusteringcoef,sizewcc" --accuracy 100 --repeat ${repeat_index}
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir es/add_global_deg_min_0.0000 --algo_opt "intersection,densification,clusteringcoef,sizewcc" --accuracy 100 --repeat ${repeat_index}
        ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helper --inputdir tihs --algo_opt "intersection,densification,sizewcc,clusteringcoef" --accuracy 100 --repeat ${repeat_index}
        cd src
        python calculation_helper.py --dataset ${data} --algorithm ns/add_global_deg_0.0000 --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
        python calculation_helper.py --dataset ${data} --algorithm ns/add_global_deg_1.0000 --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
        python calculation_helper.py --dataset ${data} --algorithm rw/rw_c_1.00_0.30 --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
        python calculation_helper.py --dataset ${data} --algorithm ff/ff_c_0.51_0.20 --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
        python calculation_helper.py --dataset ${data} --algorithm ns/add_ordered --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
        python calculation_helper.py --dataset ${data} --algorithm es/add_global_deg_min_0.0000 --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
        python calculation_helper.py --dataset ${data} --algorithm tihs --effdiam --overlapness --repeat ${repeat_index} --accuracy 100
        cd ..
        for sp in ${spset[@]}
        do
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir ns/add_global_deg_0.0000 --samplingportion ${sp} --repeat ${repeat_index}
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir ns/add_global_deg_1.0000 --samplingportion ${sp} --repeat ${repeat_index}
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir rw/rw_c_1.00_0.30 --samplingportion ${sp} --repeat ${repeat_index}
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir ff/ff_c_0.51_0.20 --samplingportion ${sp} --repeat ${repeat_index}
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir ns/add_ordered --samplingportion ${sp} --repeat ${repeat_index}
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir es/add_global_deg_min_0.0000  --samplingportion ${sp} --repeat ${repeat_index}
            ./bin/Sampling --dataname ${data} --inputpath ../dataset/ --algorithm helperdist --inputdir tihs --samplingportion ${sp} --repeat ${repeat_index}
            
            # For threads-ask-ubuntu, coauth-MAG-Geology-full, and coauth-MAG-History-full, 
            # run `src/preprocess_sv.py --repeat_str ${repeat_index} --dataname ${data}`
            # and then run `run_sv_base.m` before the next line
            
            cd src
            python calculation_helper.py --dataset ${data} --algorithm ns/add_global_deg_0.0000 --svdist --samplingportion ${sp} --repeat ${repeat_index}
            python calculation_helper.py --dataset ${data} --algorithm ns/add_global_deg_1.0000 --svdist --samplingportion ${sp} --repeat ${repeat_index}
            python calculation_helper.py --dataset ${data} --algorithm rw/rw_c_1.00_0.30 --svdist --samplingportion ${sp} --repeat ${repeat_index}
            python calculation_helper.py --dataset ${data} --algorithm ff/ff_c_0.51_0.20 --svdist --samplingportion ${sp} --repeat ${repeat_index}
            python calculation_helper.py --dataset ${data} --algorithm ns/add_ordered --svdist --samplingportion ${sp} --repeat ${repeat_index}
            python calculation_helper.py --dataset ${data} --algorithm es/add_global_deg_min_0.0000 --svdist --samplingportion ${sp} --repeat ${repeat_index}
            python calculation_helper.py --dataset ${data} --algorithm tihs --svdist --samplingportion ${sp} --repeat ${repeat_index}
            cd ..
        done
    done
done